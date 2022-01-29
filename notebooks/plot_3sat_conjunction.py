#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 19:10:03 2020

@author: blagau
"""
#PLOTS THE RESULTS
#====================================
# create fig and axes objects

# selects a data interval centered on timesB
timesFAC = [timesA[idxconj[0]], timesB[idxconj[1]], timesC[idxconj[2]]]
df_int = [dat_fac[ii][timesFAC[1] - marg : timesFAC[1] + marg] for ii in range(3)]

tbeg = dBgeo.index[0]
tend = dBgeo.index[-1]

fname_fig = 'swABC_conjunction_'+ \
        timesFAC[1].strftime('%Y%m%d_%H%M') +'.pdf'

fig_size = (8.27,11.69)
fig = plt.figure(figsize=fig_size, frameon=True)
fig.suptitle('Swarm conjunction interval: ' + tbeg.strftime('%Y-%m-%d') + \
        '   ' + tbeg.strftime('%H:%M') + ' - ' + tend.strftime('%H:%M'), \
             size='xx-large', fontweight = 'demibold', y=0.992)

# PANEL PART
# designes the panel frames 
xle, xri, ybo, yto = 0.12, 0.94, 0.03, 0.06     # plot margins
ratp = np.array([1.2, 0.6, 1.2, 0.6, 1.2, 0.6, 1., 1., 1.])         #relative hight of each panel 
hsep = np.array([0., 0., 0.07, 0., 0.07, 0., 0.07, 0.0, 0.0])        # vspace between panels 
nrp = len(ratp) 

hun = (1 - yto - ybo - hsep.sum())/ratp.sum() 
ylo = np.zeros(nrp)     # y low for each panel
yhi = np.zeros(nrp)     # y high for each panel
for ii in range(0, nrp): 
    ylo[ii] = (1 - yto) -  ratp[: ii+1].sum()*hun - hsep[: ii+1].sum()
    yhi[ii] = ylo[ii] + hun*ratp[ii]

# creates axex for each panel
ax = [0] * nrp
for ii in range(nrp):
    ax[ii] =  fig.add_axes([xle, ylo[ii], xri-xle, hun*ratp[ii]])
    
for ii in [0, 1, 2, 3, 4, 5]:    
    ax[ii].set_xlim(timesFAC[1] - marg, timesFAC[1] + marg)
    ax[ii].set_xticklabels([])  

for ii in [6, 7, 8]: 
    ax[ii].get_shared_x_axes().join(ax[6], ax[7], ax[8])


#gets the QDLat trend, sign, maximum value, and value at central FAC  
qd_trends_conj, qd_signs_conj, qd_max_conj, qd_aos_conj = \
            (np.zeros(3) for i in range(4))
for ii in range(3):
    qd_trends_conj[ii] = qd_trends_3sc[ii][idxconj[ii]]
    qd_signs_conj[ii] = qd_signs_3sc[ii][idxconj[ii]]        
    qd_max_conj[ii] = df_int[ii]['QDLat'].abs().max()
    qd_aos_conj[ii] = qd_aos_3sc[ii][idxconj[ii]]

# defines the common QDLat range for the last three panels as:
# - lower limit is the mean QDLat values of the central FACs, minus 7 degree
# - upper value as +/- maximum QDLat absolute values
# - the trend and sign is taken from SwarmA
ss = qd_signs_conj[0]
tt = qd_trends_conj[0]
qd_range = ss*np.array([abs(qd_aos_conj.mean().round()) - 7, np.ceil(qd_max_conj.max())])
if ss*tt < 0:
    qd_range = np.flip(qd_range)

# Plots title    
ax[0].set_title('\ntime differences: ' + \
   '   swA - swB = ' +  str(int((timesFAC[0] - timesFAC[1]).total_seconds())) + \
   '   swC - swB = ' +  str(int((timesFAC[2] - timesFAC[1]).total_seconds())),\
          fontsize = 'xx-large', pad = 8) 

# Plots the coresponding evolutions
for jj in range(nsc):
    df_int_sc = df_int[jj]
    # plot the fpanels with dBgeo data
    ax[jj*2].plot(dBgeo[('dBgeo',sats[jj])])
    ax[jj*2].set_ylabel('$dB_{GEOC}$ sw'+sats[jj]+'\n[nT]', linespacing=1.7)
    ax[jj*2].axvline(timesFAC[jj], ls='--', c='k')
    ax[jj*2].axvline(timesFAC[1], ls='--', c='r')
    ax[jj*2].axvline(df_int_sc['QDLat'].abs().idxmax(), ls='-', c='b')
    
    # plots the fpanels with filtered FAC data    
    ax[jj*2 + 1].plot(df_int_sc['FAC_flt_sup'], linewidth=2)
    ax[jj*2 + 1].set_ylabel('$J_{FAC}$\n[$\mu A/m^2$]', linespacing=1.7)
    ax[jj*2 + 1].axvline(timesFAC[jj], ls='--', c='k')
    ax[jj*2 + 1].axvline(timesFAC[1], ls='--', c='r')
    ax[jj*2 + 1].axvline(df_int_sc['QDLat'].abs().idxmax(), ls='-', c='b')
    ax[jj*2 + 1].axhline(0, ls='--', c='k')

    # adds QDLat, QDLon and MLT tick labels    
    locx = ax[jj*2 + 1].get_xticks()
    qdlat_ipl = np.round(np.interp(locx, mdt.date2num(df_int_sc.index), \
                            df_int_sc['QDLat']), decimals=2).astype('str')
    qdlon_ipl = np.round(np.interp(locx, mdt.date2num(df_int_sc.index), \
                            df_int_sc['QDLon']), decimals=2).astype('str')
    mlt_ipl = np.round(np.interp(locx, mdt.date2num(df_int_sc.index), \
                            df_int_sc['MLT']), decimals=1).astype('str')
    lab_fin = ['']*len(locx)
    for ix in range(len(locx)):
        lab_ini = mdt.num2date(locx[ix]).strftime('%H:%M:%S')
        lab_fin[ix] = lab_ini + '\n' +qdlat_ipl[ix] + '\n' + \
                    qdlon_ipl[ix] + '\n' + mlt_ipl[ix]
    ax[jj*2 + 1].set_xticklabels(lab_fin)
    plt.figtext(0.01, ylo[jj*2 + 1]-0.008, 'Time\nQDLat\nQDLon\nMLT', va='top')
    plt.figtext(0.96, (ylo[jj*2 + 1]+yhi[jj*2])/2, 'Swarm '+ sats[jj], va='center', \
            rotation=90, size='xx-large', fontweight = 'medium')

    # computes dBgeo as function of QDLat. For that takes the start and stop 
    # times for the quarter-orbit and applies find_jabs_midcsum to computes 
    # the evolution of time as a function of QDLat
    tbeg_qor = beg_qor_3sc[jj][idxconj[jj]]      
    tend_qor = end_qor_3sc[jj][idxconj[jj]]    
    qorb_jj = dat_fac[jj][tbeg_qor:tend_qor]    # quarter-orbit data
    tjhalf, qd_trend, qd_sign, qd_absmax, qdlon, ti_arr, jb_arr, jabs_csum, \
            qd_arr = find_jabs_midcsum(qorb_jj)
    cmp = ['X', 'Y', 'Z']
    db_arr = np.zeros((len(ti_arr),3))
    for ic in range(3):
        db_arr[:,ic] = np.interp(ti_arr, dBgeo[('dBgeo',sats[jj])].index.values.astype(float), \
                   dBgeo[('dBgeo',sats[jj],cmp[ic])].values)

    ind_qd = (np.abs(ti_arr - tjhalf)).argmin()
    # plots the last three panels with dBgeo as function of QDLat.
    ax[jj + 6].plot(qd_arr, db_arr)
    ax[jj + 6].set_ylabel('$dB_{GEOC}$ sw'+sats[jj]+'\n[nT]', linespacing=1.7)
    ax[jj + 6].axvline(qd_aos_conj[jj], ls='--', c='k')
    ax[jj + 6].axhline(0, ls='--', c='k')
    ax[jj + 6].xaxis.set_major_locator(plt.MaxNLocator(10))
    ax[jj + 6].set_xlim(qd_range)
    if jj in range(2):
        ax[jj + 6].set_xticklabels([])
        ax[jj + 6].get_shared_y_axes().join(ax[6], ax[7], ax[8])       
    plt.figtext(0.01, ylo[8]-0.008, 'QDLat', va='top')
 
fig.savefig(fname_fig)    