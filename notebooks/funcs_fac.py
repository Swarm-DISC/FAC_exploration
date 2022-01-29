#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 11:10:29 2021

@author: blagau
"""
import numpy as np
import pandas as pd
#from scipy import signal
#from scipy.interpolate import interp1d
import datetime as dtm


def split_into_sections(df, begend_arr):
    """Splits the original dataframe df in sections according to 
    time moments from begend_arr array. Returns a list of dataframes"""
    secorbs = []
    for start, end in zip(begend_arr[0:-1], begend_arr[1:]):
        secorb = df[start:end]
        secorbs.append(secorb)
    return secorbs

def normvec(v):
    # computes the unit-vectors of a vector time series
    return np.divide(v,np.linalg.norm(v,axis=-1).reshape(-1,1))


def rotvecax(v, ax, ang):
    # Rotates vector v by angle ang around vector ax 
    # Uses Rodrigues' formula when v is normal to ax
    sa, ca = np.sin(np.deg2rad(ang)), np.cos(np.deg2rad(ang))
    return v*ca[...,np.newaxis] + np.cross(ax, v)*sa[...,np.newaxis]


def sign_ang(V, N, R):
    # returns the signed angle between vectors V and N, perpendicular 
    # to R; positive sign corresponds to right hand rule along R
    VxN = np.cross(V, N)
    pm = np.sign(np.sum(R*VxN, axis=-1))
    return np.degrees(np.arcsin(pm*np.linalg.norm(VxN, axis=-1)))       
    

def R_B_dB_in_GEOC(Rsph, Bnec, dBnec):
    # transforms magnetic field and magnetic field perturbation 
    # from NEC to GEO. Transform position from spherical 
    # to GEO cartesian
    latsc = np.deg2rad(Rsph[:,0])
    lonsc = np.deg2rad(Rsph[:,1])  
    radsc = 0.001*Rsph[:,2]
    # prepares conversion to global cartesian frame
    clt,slt = np.cos(latsc.flat),np.sin(latsc.flat)
    cln,sln = np.cos(lonsc.flat),np.sin(lonsc.flat)
    north = np.stack((-slt*cln,-slt*sln,clt),axis=-1)
    east = np.stack((-sln,cln,np.zeros(cln.shape)),axis=-1)
    center = np.stack((-clt*cln,-clt*sln,-slt),axis=-1)
    # stores cartesian position vectors in position data matrix R
    R = -radsc[...,None]*center
    # stores magnetic data in GEOC (same frame as for R)
    Bgeo = np.matmul(np.stack((north,east,center),axis=-1),
                        Bnec[...,None]).reshape(Bnec.shape)
    dBgeo = np.matmul(np.stack((north,east,center),axis=-1),
                        dBnec[...,None]).reshape(dBnec.shape)
    return R, Bgeo, dBgeo


def GapsAsNaN(df_ini, ind_gaps):
    # Replaces with NaN the zero values in L1b 
    # low-resolution magnetic field files.
    df_out = df_ini.copy()
    df_out['B_NEC'][ind_gaps] = [np.full(3,np.NAN)]*len(ind_gaps)
    return df_out, df_ini.index[ind_gaps]


def singleJfac(t, R, B, dB, alpha=None, N2d=None, \
               N3d=None, tincl=None, er_db=0.5):
    # Computes single-satellite FAC density 
    
    # Constructs the differences & values at mid-intervals
    dt = t[1:].values - t[:-1].values
    tmid = t[:-1].values + dt*0.5
    Bmid = 0.5*(B[1:,:] + B[:-1,:])           
    Rmid = 0.5*(R[1:,:] + R[:-1,:])
    diff_dB = dB[1:,:] - dB[:-1,:]    
    V3d = R[1:,:] - R[:-1,:]  
    Vorb = np.sqrt(np.sum(V3d*V3d, axis=-1))      
    # Defines important unit vectors
    eV3d, eBmid, eRmid = normvec(V3d), normvec(Bmid), normvec(Rmid)
    eV2d = normvec(np.cross(eRmid, np.cross(eV3d, eRmid)))    
    # Angle between B and R
    cos_b_r = np.sum(eBmid*eRmid, axis=-1)
    bad_ang = np.abs(cos_b_r) < np.cos(np.deg2rad(60))
 
    # incl is the array of FAC incliation wrt Vsat (in tangential plane)    
    if N3d is not None:
        eN3d = normvec(N3d)
        eN2d = normvec(eN3d - np.sum(eN3d*eRmid,axis=-1).reshape(-1,1)*eRmid)
        incl = sign_ang(eV2d, eN2d, eRmid)
    elif alpha is not None:
        incl = alpha if isinstance(alpha, np.ndarray) else \
                                     np.full(len(tmid), alpha)        
    elif N2d is not None:
        eN2d = normvec(np.cross(eRmid, np.cross(N2d, eRmid)))
        incl = sign_ang(eV2d, eN2d, eRmid)
    else:
        incl = np.zeros(len(tmid))

    # considers the validity interval of FAC inclination 
    if tincl is not None:
        ind_incl = np.where((tmid >= tincl[0]) & (tmid <= tincl[1]))[0]
        incl[0:ind_incl[0]] = incl[ind_incl[0]]
        incl[ind_incl[-1]:] = incl[ind_incl[-1]]

    # working in the tangential plane
    eNtang = normvec(rotvecax(eV2d, eRmid, incl))
    eEtang = normvec(np.cross(eNtang, eRmid))
    diff_dB_Etang = np.sum(diff_dB*eEtang, axis=-1)
    Dplane = np.sum(eNtang*eV2d, axis=-1)
    j_rad= - diff_dB_Etang/Dplane/Vorb/(4*np.pi*1e-7)*1.e-6
    j_rad_er= np.abs(er_db/Dplane/Vorb/(4*np.pi*1e-7)*1.e-6)   
    
    # FAC density and error
    j_b = j_rad/cos_b_r
    j_b_er = np.abs(j_rad_er/cos_b_r)    
    j_b[bad_ang] = np.nan
    j_b_er[bad_ang] = np.nan    
    
    return tmid, Rmid, j_b, j_rad, j_b_er, j_rad_er, incl, np.arccos(cos_b_r)*180./np.pi


def find_ao_margins(df, fac_qnt = 'FAC_flt_sup', rez_qd = 100):
    """For a quarter-orbit section, approximate the 
    margins of auroral oval mar"""
    qd = df['QDLat'].values
    jb = df[fac_qnt].values
    ti = df['QDLat'].index.values.astype(float)
    qd_trend = (qd[-1] - qd[0])/abs(qd[0] - qd[-1])  # >0 if QDLat inreases
    qd_sign = qd[0]/abs(qd[0])  
    dqmin = np.round(np.amin(qd), decimals = 2)
    dqmax = np.round(np.amax(qd), decimals = 2) 
    nr = round((dqmax - dqmin)*rez_qd + 1)      
    qd_arr = np.linspace(dqmin, dqmax, nr)       # new QDLat array 
    if qd_trend > 0:
        ti_arr = np.interp(qd_arr, qd, ti)          # time of new QD points
        jb_arr = np.interp(qd_arr, qd, jb)          # FAC at new QD points
    else:
        ti_arr = np.interp(qd_arr, qd[::-1], ti[::-1])
        jb_arr = np.interp(qd_arr, qd[::-1], jb[::-1])
    
    if np.sum(np.abs(jb_arr)) > 0:
        jabs_csum = np.cumsum(np.abs(jb_arr))/np.sum(np.abs(jb_arr))
        idx_j1q = (np.abs(jabs_csum - 0.25)).argmin()
        idx_j2q = (np.abs(jabs_csum - 0.5)).argmin()
        idx_j3q = (np.abs(jabs_csum - 0.75)).argmin()
        if qd_trend < 0:
            idx_j1q, idx_j3q = idx_j3q, idx_j1q
        idx_beg = np.max(np.array([0, idx_j1q - (idx_j3q - idx_j1q)]))
        idx_end = np.min(np.array([idx_j3q + (idx_j3q - idx_j1q), nr-1]))
    
        t64 = pd.DatetimeIndex(ti_arr[[idx_beg, idx_j2q, idx_end]]) 
        return t64[0], t64[1], t64[2], qd_trend, \
                qd_sign, qd_arr[idx_j2q], ti_arr, jb_arr, jabs_csum, qd_arr  
    
    else:
        return np.nan, np.nan, np.nan, qd_trend, qd_sign, np.nan, \
                ti_arr, jb_arr, np.nan, qd_arr
    
