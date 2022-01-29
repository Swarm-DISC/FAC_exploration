tq_beg = min([qorbs_dB[0][jj].index[0], qorbs_dB[1][jj].index[0]])
tq_end = max([qorbs_dB[0][jj].index[-1], qorbs_dB[1][jj].index[-1]])

iref = iref_arr[jj]
isec = (iref +1) % 2
sc_ref = sats[iref]
sc_sec = sats[isec]

# quarter orbit time and data for the second and reference s/c
tsec = qorbs_mva[isec][jj].index
dBsec = qorbs_mva[isec][jj]['dB_max'].values

tref = qorbs_mva[iref][jj].index
dBref = qorbs_mva[iref][jj]['dB_max'].values
    
str_trange = tq_beg.isoformat().replace('-','')[:8]+'_'+ \
     tq_beg.isoformat().replace(':','')[11:17] + \
        '_'+tq_end.isoformat().replace(':','')[11:17]

fname_fig = 'QI_sw'+sc_ref+sc_sec+'_'+ str_trange +'.pdf'

fig_size = (8.27,11.69)
fig = plt.figure(figsize=fig_size, frameon=True)
fig.suptitle('Quality indices for sw'+sats[0]+'/sw'+sats[1]+' orbit '+ \
             str(int(orbs[0][0]) + int(jj/4))+ '/' +str(int(orbs[1][0]) + int(jj/4))+  \
             '  quadrant ' + str(jj - 4*int(jj/4) + 1) +
             '\n', size='xx-large', y=0.995)

# MVA and cc results
all_int = 'Time interval: ' + tq_beg.isoformat()[:10]+'  ' +\
              tq_beg.isoformat()[11:19] + ' - '+tq_end.isoformat()[11:19]

all_cc = 'Correlation analysis (ref. sat.,  coeff,  time lag [s]):    sw' + \
            sats[iref] + '     ' + str(np.round(cc_ls[jj], decimals=3)) + \
        '      '+ str(opt_lag_ls[jj]) 

cap_mva = 'MVA results (sat,  interval,  $\lambda_{max}/\lambda_{min}$,  ' +\
            'angle V N [deg.])'
sw0 = '      sw'+sats[0]+':   '+ tbeg_mva[0][jj].isoformat()[11:19] + \
        ' - '+tbeg_mva[0][jj].isoformat()[11:19] + '      ' +\
        str(np.round(lbd_max[0, jj]/lbd_min[0, jj], decimals=1)) + '     ' + \
        str(np.round(ang_vn[0, jj], decimals=1))
sw1 = '      sw'+sats[1]+':   '+ tbeg_mva[1][jj].isoformat()[11:19] + \
        ' - '+tbeg_mva[1][jj].isoformat()[11:19] + '      ' +\
        str(np.round(lbd_max[1, jj]/lbd_min[1, jj], decimals=1)) + '     ' + \
        str(np.round(ang_vn[1, jj], decimals=1)) 


plt.figtext(0.25, 0.955, all_int, size='large')

plt.figtext(0.1, 0.93, cap_mva, size='large')
plt.figtext(0.1, 0.91, sw0, size='large')
plt.figtext(0.1, 0.89, sw1, size='large')
plt.figtext(0.1, 0.865, all_cc, size='large')

# PANEL PART
# designes the panel frames 
xle, xri, ybo, yto = 0.125, 0.92, 0.095, 0.15     # plot margins
ratp = np.array([1, 1, 1, 1, 1, 1, 0.6 ]) # relative hight of each panel 
hsep = 0.005                                    # vspace between panels 
nrp = len(ratp) 

hun = (1-yto - ybo - (nrp-1)*hsep)/ratp.sum() 
yle = np.zeros(nrp)     # y left for each panel
yri = np.zeros(nrp) 
for ii in range(nrp): 
    yle[ii] = (1 - yto) -  ratp[: ii+1].sum()*hun - ii*hsep
    yri[ii] = yle[ii] + hun*ratp[ii]

# creates axex for each panel
ax = [0] * nrp
for ii in range(nrp):
    ax[ii] =  fig.add_axes([xle, yle[ii], xri-xle, hun*ratp[ii]])
    ax[ii].set_xlim(tq_beg, tq_end)
    
for ii in range(nrp -1):
    ax[ii].set_xticklabels([])
    
#Plot time-series quantities    

ax[0].plot(qorbs_dB[0][jj]['dBnec'])
ax[0].legend(['dB_N', 'dB_E', 'dB_C' ], loc = (0.95, 0.1), handlelength=1)
ax[0].set_ylabel('$dB_{NEC}$ sw'+sats[0]+'\n[nT]', linespacing=1.7)
ax[0].axvline(tbeg_mva[0][jj], ls=':', c='k', lw=0.7)
ax[0].axvline(tend_mva[0][jj], ls=':', c='k', lw=0.7)

ax[1].plot(qorbs_dB[1][jj]['dBnec'])
ax[1].legend(['dB_N', 'dB_E', 'dB_C' ], loc = (0.95, 0.1), handlelength=1)
ax[1].set_ylabel('$dB_{NEC}$ sw'+sats[1]+'\n[nT]', linespacing=1.7)
ax[1].get_shared_y_axes().join(ax[0], ax[1])
ax[1].axvline(tbeg_mva[1][jj], ls=':', c='k', lw=0.7)
ax[1].axvline(tend_mva[1][jj], ls=':', c='k', lw=0.7)

ax[2].plot(qorbs_mva[0][jj][['dB_max','dB_min','dB_B']])
ax[2].legend(['dB_max', 'dB_min', 'dB_B' ], loc = (0.95, 0.1), handlelength=1)
ax[2].set_ylabel('$dB_{MVA}$ sw'+sats[0]+'\n[nT]', linespacing=1.7)

ax[3].plot(qorbs_mva[1][jj][['dB_max','dB_min','dB_B']])
ax[3].legend(['dB_max', 'dB_min', 'dB_B' ], loc = (0.95, 0.1), handlelength=1)
ax[3].set_ylabel('$dB_{MVA}$ sw'+sats[1]+'\n[nT]', linespacing=1.7)
ax[3].get_shared_y_axes().join(ax[2], ax[3])

ax[4].plot(tsec, dBsec)
qorbs_ref_cc[jj].plot(ax=ax[4])
ax[4].legend(['dB_max_sw'+sats[isec], 'dB_max_sw'+sats[iref]], \
                         loc = (0.9, 0.1), handlelength=1)
ax[4].set_ylabel('$dB_{MVA}$ \n[nT]', linespacing=1.7)

qorbs_fac[0][jj]['FAC_flt'].plot(ax=ax[5])
qorbs_fac[1][jj]['FAC_flt'].plot(ax=ax[5])
ax[5].legend(['Jfac_sw'+sats[0], 'Jfac_sw'+sats[1]], \
                         loc = (0.94, 0.1), handlelength=1)
ax[5].set_ylabel('$J_{FAC}$ \n[nT]', linespacing=1.7)

ax[6].plot(qorbs_mva[0][jj]['ang_v_n'])
ax[6].plot(qorbs_mva[1][jj]['ang_v_n'])
max_ang= max([qorbs_mva[0][jj]['ang_v_n'].max(),qorbs_mva[1][jj]['ang_v_n'].max()])
min_ang= min([qorbs_mva[0][jj]['ang_v_n'].min(),qorbs_mva[1][jj]['ang_v_n'].min()])
off_ang = 0.2*(max_ang - min_ang)
ax[6].set_ylim([min_ang - off_ang, max_ang + off_ang])
ax[6].legend(['ang_sw'+sats[0], 'ang_sw'+sats[1]], loc = (0.95, 0.1), handlelength=1)
ax[6].set_ylabel('ang_V_N \n[deg]', linespacing=1.7)
ax[6].axvline(tbeg_mva[0][jj], ls=':', c='b', lw=0.7)
ax[6].axvline(tend_mva[0][jj], ls=':', c='b', lw=0.7)
ax[6].axvline(tbeg_mva[1][jj], ls=':', c='tab:orange', lw=0.7)
ax[6].axvline(tend_mva[1][jj], ls=':', c='tab:orange', lw=0.7)
int_min = 2
dif_min = tq_end.floor('min') - tq_beg.ceil('min')
nr_ticks = dif_min/pd.to_timedelta(1, unit='m')//int_min +1
pos_ticks = tq_beg.ceil('min') + np.arange(nr_ticks)*dtm.timedelta(minutes=int_min)
ax[6].set_xticks(pos_ticks)


# Ephemerides
latc = qorbs_Bnec[0][jj]['Latitude'].values
lonc = qorbs_Bnec[0][jj]['Longitude'].values

locx = ax[nrp-1].get_xticks()
latc_ipl = np.round(np.interp(locx, mdt.date2num(qorbs_mva[0][jj].index), \
                            latc), decimals=2).astype('str')
lonc_ipl = np.round(np.interp(locx, mdt.date2num(qorbs_mva[0][jj].index), \
                            lonc), decimals=2).astype('str')
qdlat_ipl = np.round(np.interp(locx, mdt.date2num(qorbs_mva[0][jj].index), \
                            qorbs_Bnec[0][jj]['QDLat']), decimals=2).astype('str')
qdlon_ipl = np.round(np.interp(locx, mdt.date2num(qorbs_mva[0][jj].index), \
                            qorbs_Bnec[0][jj]['QDLon']), decimals=2).astype('str')
mlt_ipl = np.round(np.interp(locx, mdt.date2num(qorbs_mva[0][jj].index), \
                            qorbs_Bnec[0][jj]['MLT']), decimals=1).astype('str')

lab_fin = ['']*len(locx)
for ii in range(len(locx)):
    lab_ini = mdt.num2date(locx[ii]).strftime('%H:%M:%S')
    lab_fin[ii] = lab_ini + '\n' +latc_ipl[ii] + '\n' + lonc_ipl[ii] + \
    '\n'+ qdlat_ipl[ii] + '\n' +qdlon_ipl[ii] + '\n' + mlt_ipl[ii]
    
ax[nrp-1].set_xticklabels(lab_fin)
plt.figtext(0.01, 0.01, 'Time\nLat\nLon\nQLat\nQLon\nMLT')


plt.show()
fig.savefig(fname_fig)