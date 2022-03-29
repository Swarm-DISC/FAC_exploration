#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 11:09:38 2021

@author: blagau
"""

#SAVES INPUT AND OUTPUT TO ASCII FILES
#====================================
str_trange = dtime_beg.replace('-','')[:8]+'_'+ \
     dtime_beg.replace(':','')[11:17] + '_'+dtime_end.replace(':','')[11:17]

fname_out = 'FAC_'+rez+'_sw'+sat[0]+'_atmuso_'+ str_trange +'.dat'
fname_fig = 'FAC_'+rez+'_sw'+sat[0]+'_atmuso_'+ str_trange +'.pdf'

bad_ang = np.less(np.abs(np.cos(ang_BR*np.pi/180.)), np.cos(np.deg2rad(60)))

if rez == 'LR':
    info_points = '\nres.: LR     missing points: '+str(len(ind_gaps)) + '      ' +\
            'low-latitude points: ' +str(len(np.where(bad_ang)[0]))
else:
    info_points = '\nres.: HR     low-latitude points: '  +str(len(np.where(bad_ang)[0]))                
add_text = 'Time interval:  ' + dtime_beg[:10]+'  ' + \
            dtime_beg[11:19] + ' - '+dtime_end[11:19] + info_points
           
if tincl is not None:
    add_text = add_text + '\nInclination interval:  ' + \
        np.datetime_as_string(tincl[0])[11:19] + ' - ' + \
        np.datetime_as_string(tincl[1])[11:19]
    
# exports the results
out_df = j_df.copy()   
text_trj = '# time [YYYY-mm-ddTHH:MM:SS.f],\t  Rmid_X [km],\t  Rmid_Y,\t\
  Rmid_Z,\t  Jfac [microA/m^2],\t  IRC,\t  errJfac,\t  errJirc,\t  '
    
if use_filter:
    out_df = pd.concat([j_df[['Rmid X','Rmid_Y','Rmid Z','FAC','IRC','FAC_er','IRC_er']],\
         jflt_df[['FAC_flt','IRC_flt','FAC_flt_er','IRC_flt_er','ang_BR','incl']]], axis=1)
    text_trj = text_trj + 'Jfac_flt,\t  Jirc_flt,\t  errJfac_flt,\t  errJirc_flt,\t  '         
text_header = text_trj + 'ang_BR [deg],\t  incl \n'

with open(fname_out, 'w') as file:
    file.write('# Swarm sat: ' + str(sat[0]) + '\n')
    file.write('# Model: ' + Bmodel + '\n')
    if rez == 'LR':
        file.write('# Number of missing data points in L1b file: ' + str(len(ind_gaps)) + '\n')
        if len(ind_gaps) > 0:
            file.write('#' + '; '.join(timegaps.strftime('%H:%M:%S').values) + '\n')
    file.write(text_header)
out_df.to_csv(fname_out, mode='a', sep=",", na_rep = 'NaN', date_format='%Y-%m-%dT%H:%M:%S.%f',\
              float_format='%15.4f', header=False)


#PLOTS THE RESULTS
#====================================
# creates fig and axes objects
fig_size = (8.27,11.69)
fig = plt.figure(figsize=fig_size, frameon=True)
fig.suptitle('FAC density estimate with single-satellite method\n', size='xx-large',\
             weight = 'bold', y=0.99)
plt.figtext(0.21, 0.96, add_text, fontsize = 'x-large', \
            va='top', ha = 'left')

xle, xri, ybo, yto = 0.125, 0.92, 0.095, 0.105     # plot margins
ratp = np.array([1, 1, 1, 1, 0.5, 0.5])          #relative hight of each panel 
hsep = 0.005                                    # vspace between panels 
nrp = len(ratp) 

hun = (1-yto - ybo - (nrp-1)*hsep)/ratp.sum() 
yle = np.zeros(nrp)     # y left for each panel
yri = np.zeros(nrp) 
for ii in range(nrp): 
    yle[ii] = (1 - yto) -  ratp[: ii+1].sum()*hun - ii*hsep
    yri[ii] = yle[ii] + hun*ratp[ii]

# creates axes for each panel
ax = [0] * nrp
for ii in range(nrp):
    ax[ii] =  fig.add_axes([xle, yle[ii], xri-xle, hun*ratp[ii]])
    ax[ii].set_xlim(pd.Timestamp(dtime_beg), pd.Timestamp(dtime_end))
  
for ii in range(nrp -1):
    ax[ii].set_xticklabels([])
    
#Plots time-series quantities    
ax[0].plot(dB_df)
ax[0].set_ylabel('$dB_{GEO}$ sw'+sat[0]+'\n[nT]', linespacing=1.7)
ax[0].legend(['dB_X', 'dB_Y', 'dB_Z' ], loc='upper right')

ax[1].plot(j_df['FAC'], label='FAC')
if use_filter:
    ax[1].plot(jflt_df['FAC_flt'], label='FAC_flt')
ax[1].get_shared_x_axes().join(ax[0], ax[1], ax[2], ax[3], ax[4], ax[5])
ax[1].set_ylabel(r'$J_{FAC}$'+'\n'+r'$[\mu A/m^2]$', linespacing=1.7)
ax[1].legend(loc='upper right')

ax[2].plot(j_df['IRC'], label='IRC')
if use_filter:
    ax[2].plot(jflt_df['IRC_flt'], label='IRC_flt')
ax[2].get_shared_x_axes().join(ax[0], ax[1], ax[2], ax[3], ax[4], ax[5])
ax[2].set_ylabel(r'$J_{IRC}$'+'\n'+r'$[\mu A/m^2]$', linespacing=1.7)
ax[2].legend(loc='upper right') 

ax[3].plot(j_df['FAC'], label='FAC')
ax[3].plot(FAC_L2['FAC_L2'], label='FAC_L2')
ax[3].get_shared_x_axes().join(ax[0], ax[1], ax[2], ax[3], ax[4], ax[5])
ax[3].set_ylabel(r'$J_{FAC}$'+'\n'+r'$[\mu A/m^2]$', linespacing=1.7)
ax[3].legend(loc='upper right') 

ax[4].plot(j_df['ang_BR'])
ax[4].get_shared_x_axes().join(ax[0], ax[1], ax[2], ax[3], ax[4], ax[5])
ax[4].set_ylabel('ang. BR'+'\n'+r'$[deg]$', linespacing=1.7)
ax[4].legend() 

ax[5].plot(j_df['incl'])
ax[5].get_shared_x_axes().join(ax[0], ax[1], ax[2], ax[3], ax[4], ax[5])
ax[5].set_ylabel('incl.'+'\n'+r'$[deg]$', linespacing=1.7)
ax[5].legend() 

if tincl is not None:
    for ii in range(nrp):
        ax[ii].axvline(tincl[0], ls='--', c='k')
        ax[ii].axvline(tincl[1], ls='--', c='k')

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
plt.figtext(0.01, 0.01, 'Time\nLat\nLon\nQLat\nQLon\nMLT')


plt.show()

fig.savefig(fname_fig)