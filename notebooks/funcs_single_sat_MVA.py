#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 26 11:10:29 2021

@author: blagau
"""
import numpy as np
#import pandas as pd
#from scipy import signal
#from scipy.interpolate import interp1d
#import datetime as dtm


def normvec(v):
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
    latsc = np.deg2rad(Rsph[:,0])
    lonsc = np.deg2rad(Rsph[:,1])  
    radsc = 0.001*Rsph[:,2]
    # prepare conversion to global cartesian frame
    clt,slt = np.cos(latsc.flat),np.sin(latsc.flat)
    cln,sln = np.cos(lonsc.flat),np.sin(lonsc.flat)
    north = np.stack((-slt*cln,-slt*sln,clt),axis=-1)
    east = np.stack((-sln,cln,np.zeros(cln.shape)),axis=-1)
    center = np.stack((-clt*cln,-clt*sln,-slt),axis=-1)
    # store cartesian position vectors in position data matrix R
    R = -radsc[...,None]*center
    # store magnetic data in GEOC (same frame as for R)
    Bgeo = np.matmul(np.stack((north,east,center),axis=-1),
                        Bnec[...,None]).reshape(Bnec.shape)
    dBgeo = np.matmul(np.stack((north,east,center),axis=-1),
                        dBnec[...,None]).reshape(dBnec.shape)
    return R, Bgeo, dBgeo

# Replaces with NaN the zero values in L1b low-resolution magnetic field files.
def GapsAsNaN(df_ini, ind_gaps):
    df_out = df_ini.copy()
    df_out['B_NEC'][ind_gaps] = [np.full(3,np.NAN)]*len(ind_gaps)
    return df_out, df_ini.index[ind_gaps]

# function that computes single-satellite FAC density
def singleJfac(t, R, B, dB, alpha=None, N2d=None, \
               N3d=None, tincl=None, er_db=0.5):
    # Construct the differences & values at mid-intervals
    dt = t[1:].values - t[:-1].values
    tmid = t[:-1].values + dt*0.5
    Bmid = 0.5*(B[1:,:] + B[:-1,:])           
    Rmid = 0.5*(R[1:,:] + R[:-1,:])
    diff_dB = dB[1:,:] - dB[:-1,:]    
    V3d = R[1:,:] - R[:-1,:]  
    Vorb = np.sqrt(np.sum(V3d*V3d, axis=-1))      
    # Define important unit vectors
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

    # consider the validity interval of FAC inclination 
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

