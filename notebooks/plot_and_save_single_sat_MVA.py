#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 11:09:38 2021

@author: blagau
"""

#SAVE INPUT AND OUTPUT TO ASCII FILES
#====================================
trange = [tmva_int[0].astype(dtm.datetime) - tmarg, \
                          tmva_int[1].astype(dtm.datetime) + tmarg]

str_trange = trange[0].strftime("%Y-%m-%d   %H:%M:%S") +\
           ' - ' + trange[1].strftime("%H:%M:%S")

str_tmva = (tmva_int[0].astype(dtm.datetime)).strftime("%Y-%m-%d   %H:%M:%S") +\
           ' - ' + (tmva_int[1].astype(dtm.datetime)).strftime("%H:%M:%S")

tmva_mid = tmva_int[0] + (tmva_int[1] - tmva_int[0])/2.

str_fname = (tmva_int[0].astype(dtm.datetime)).strftime("%Y%m%d_%H%M%S") +\
           '_' + (tmva_int[1].astype(dtm.datetime)).strftime("%H%M%S")


indfile = np.where((ti >= trange[0]) & (ti <= trange[1]))[0]

fname_out = 'MVA_sw'+sat[0] + '_'+ str_fname +'.dat'
fname_fig = 'MVA_sw'+sat[0] + '_'+ str_fname +'.eps'

mva_results = 'MVA interval:  '+ \
     tmva_int[0].astype(dtm.datetime).strftime("%H:%M:%S") + ' - ' + \
     tmva_int[1].astype(dtm.datetime).strftime("%H:%M:%S") +'\n' + \
    r'$\mathbf{N_{GEO}}$' +':'+ '  ['+ '{:.4f},  '.format(mindir[0]) +\
            '{:.4f},  '.format(mindir[1]) + '{:.4f}]'.format(mindir[2]) +'\n' + \
    'eigenvalues ratio: ' + '{:.2f}'.format(eigval[2]/eigval[1]) +'\n' + \
    'FAC inclination (positive from '+r'$\mathbf{V_{sat}}$'+ ' to '+ \
    r'$\mathbf{N_{tang}}$' + ' along '+ r'$\mathbf{R_{sat}}$' + '):  ' + \
    '{:.1f}'.format(np.min(ang[indok[:-1]])) + r'$^{\degree}$' +\
    ' ' +r'$\mathrm{\div}$'+ '   ' + '{:.1f}'.format(np.max(ang[indok[:-1]])) + r'$^{\degree}$'       
         
    
# exports the results
with open(fname_out, 'w') as file:
    file.write('# Swarm sat: ' + str(sat[0]) + '\n')
    file.write('# Model: ' + Bmodel + '\n')
    file.write('# time interval: ' + str_trange + '\n')   
    file.write('#########  MVA results  #########\n')
    file.write('# MVA interval: ' + str_tmva + '\n')    
    file.write('# B_unit:  '+ format(eigval[0],'.2f')+'  ['+ format(B_unit[0],'.4f')\
        + ', '+format(B_unit[1], '.4f') +', '+format(B_unit[2], '.4f') +'] \n')    
    file.write('# minvar:  '+ format(eigval[1],'.2f')+'  ['+ format(mindir[0],'.4f')\
        + ', '+format(mindir[1], '.4f') +', '+format(mindir[2], '.4f') +'] \n')  
    file.write('# maxvar:  '+ format(eigval[2],'.2f')+'  ['+ format(maxdir[0],'.4f')\
        + ', '+format(maxdir[1], '.4f') +', '+format(maxdir[2], '.4f') +'] \n') 
    file.write('# eigenvalues ratio: ' + format(eigval[2]/eigval[1], '.1f') +'\n')    
    file.write('# FAC inclination wrt Vsat (tangential plane):  ' + \
        format(np.min(ang[indok[:-1]]), '.1f') + ' -  ' + \
        format(np.max(ang[indok[:-1]]), '.1f') + '  deg. \n') 
    file.write('##############################\n') 
    file.write('# time [YYYY-mm-ddTHH:MM:SS],\tdB_mva_{B, minvar, maxvar} [nT] \n')

dBmva_df.iloc[indfile].to_csv(fname_out, mode='a', sep=",", date_format='%Y-%m-%dT%H:%M:%S',\
              float_format='%15.3f', header=False)


#PLOTS THE RESULTS
#====================================
# creates fig and axes objects
fig_size = (8.27,11.69)
fig = plt.figure(figsize=fig_size, frameon=True)
fig.suptitle('MVA results on Swarm'+sat[0] + '  ' +\
    tmva_mid.astype(dtm.datetime).strftime("%Y-%m-%d  %H:%M:%S"), \
    size='xx-large', weight = 'bold', y=0.99)

plt.figtext(0.04, 0.95, mva_results, fontsize = 'x-large', \
            va='top', ha = 'left')

xle, xri, ybo, yto = 0.125, 0.94, 0.45, 0.15     # plot margins
ratp = np.array([1, 1, 1, 0.5]) #relative hight of each panel 
hsep = 0.005                                    # vspace between panels 
nrp = len(ratp) 

hun = (1-yto - ybo - (nrp-1)*hsep)/ratp.sum() 
yle = np.zeros(nrp)     # y left for each panel
yri = np.zeros(nrp) 
for ii in range(nrp): 
    yle[ii] = (1 - yto) -  ratp[: ii+1].sum()*hun - ii*hsep
    yri[ii] = yle[ii] + hun*ratp[ii]

# creates axes for each panel
ax = [0] * (nrp+1)
for ii in range(nrp):
    ax[ii] =  fig.add_axes([xle, yle[ii], xri-xle, hun*ratp[ii]])
    ax[ii].set_xlim(pd.Timestamp(trange[0]), pd.Timestamp(trange[1]))

for ii in range(nrp -1):
    ax[ii].set_xticklabels([])
    
#Plots time-series quantities    
ax[0].plot(dBnec_df)
ax[0].set_ylabel('$dB_{GEO}$ sw'+sat[0]+'\n[nT]', linespacing=1.7)
ax[0].axvline(tmva_int[0], ls='--', c='k')
ax[0].axvline(tmva_int[1], ls='--', c='k')
ax[0].legend(['dB_N', 'dB_E', 'dB_C' ], loc = (0.95, 0.1), handlelength=1)
ax[0].axhline(0, ls='--', c='k')

ax[1].plot(dBmva_df)
ax[1].set_ylabel('$dB_{MVA}$ sw'+sat[0]+'\n[nT]', linespacing=1.7)
ax[1].axvline(tmva_int[0], ls='--', c='k')
ax[1].axvline(tmva_int[1], ls='--', c='k')
ax[1].legend(['dB_B', 'dB_min', 'dB_max' ], loc = (0.93, 0.1), handlelength=1)
ax[1].axhline(0, ls='--', c='k')

if use_filter:
    ax[2].plot(jb_flt_df, label='FAC')
    ax[2].plot(jb_inc_flt_df, label='FAC_inc')
else:
    ax[2].plot(jb_df, label='FAC')
    ax[2].plot(jb_inc_df, label='FAC_inc')    
ax[2].axhline(0, ls='--', c='k')
ax[2].axvline(tmva_int[0], ls='--', c='k')
ax[2].axvline(tmva_int[1], ls='--', c='k')
ax[2].set_ylabel(r'$J_{FAC}$'+'\n'+r'$[\mu A/m^2]$', linespacing=1.7)
ax[2].legend(loc = (0.93, 0.6), handlelength=1) 

ax[3].plot(ang_df)
ax[3].axvline(tmva_int[0], ls='--', c='k')
ax[3].axvline(tmva_int[1], ls='--', c='k')
ax[3].set_ylabel(r'$ang_NV$'+'\n'+r'$[deg]$', linespacing=1.7)

ax[4] =  fig.add_axes([xle, 0.03, xri-xle, 0.30])
ax[4].set_title('Hodogram of '+r'$dB_{minvar}$'+ ' vs. '+\
    r'$dB_{maxvar}$' , fontsize = 'xx-large', pad = 15)
ax[4].plot(dBmva_df.values[indfile,2], dBmva_df.values[indfile,1], label='plot_range')
ax[4].plot(dBmva_df.values[indok,2], dBmva_df.values[indok,1], label='MVA range')
ax[4].plot(dBmva_df.values[indok[0],2], dBmva_df.values[indok[0],1], \
  label='start', color='green', marker='o', linewidth=2, markersize=8)
ax[4].plot(dBmva_df.values[indok[-1],2], dBmva_df.values[indok[-1],1], \
  label='stop', color='red', marker='o', linewidth=2, markersize=8)
ax[4].set_aspect('equal', adjustable='box')
ax[4].set_xlabel(r'$dB_{maxvar}$'+'\n'+r'$[nT]$', linespacing=1.7)
ax[4].set_ylabel(r'$dB_{minvar}$'+'\n'+r'$[nT]$', linespacing=1.7)
ax[4].legend(loc = (0.9, 0.9), handlelength=1) 


latc = Rsph[:,0]
lonc = Rsph[:,1]

locx = ax[nrp-1].get_xticks()
latc_ipl = np.round(np.interp(locx, mdt.date2num(dat_df.index), \
                            latc), decimals=2).astype('str')
lonc_ipl = np.round(np.interp(locx, mdt.date2num(dat_df.index), \
                            lonc), decimals=2).astype('str')
qdlat_ipl = np.round(np.interp(locx, mdt.date2num(dat_df.index), \
                            dat_df['QDLat']), decimals=2).astype('str')
qdlon_ipl = np.round(np.interp(locx, mdt.date2num(dat_df.index), \
                            dat_df['QDLon']), decimals=2).astype('str')
mlt_ipl = np.round(np.interp(locx, mdt.date2num(dat_df.index), \
                            dat_df['MLT']), decimals=1).astype('str')

lab_fin = ['']*len(locx)
for ii in range(len(locx)):
    lab_ini = mdt.num2date(locx[ii]).strftime('%H:%M:%S')
    lab_fin[ii] = lab_ini + '\n' +latc_ipl[ii] + '\n' + lonc_ipl[ii]+ \
    '\n'+ qdlat_ipl[ii] + '\n' +qdlon_ipl[ii] + '\n' + mlt_ipl[ii]
    
ax[nrp-1].set_xticklabels(lab_fin)
plt.figtext(0.01, 0.368, 'Time\nLat\nLon\nQLat\nQLon\nMLT')


plt.show()

fig.savefig(fname_fig)