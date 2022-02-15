#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 12:05:43 2020

@author: blagau
"""

#SAVES INPUT AND OUTPUT TO ASCII FILES
#====================================
str_trange = dtime_beg.replace('-','')[:8]+'_'+ \
     dtime_beg.replace(':','')[11:17] + '_'+dtime_end.replace(':','')[11:17]

ss = sats[0]+sats[1]
fname_in = 'input_dual_sat_atmuso_' + str_trange +'.dat'
fname_out = 'FAC_sw'+ss+'_atmuso_'+ str_trange +'.dat'
# fname_FlagsB = 'FlagsB_nonzero_'+ str_trange +'.dat'
fname_fig = 'FAC_sw'+ss+'_atmuso_'+ str_trange +'.pdf'
# export the results
namecol =  ['Rcx','Rcy','Rcz','ux ','uy ','uz ', 'angBN ', \
            'log_CN ', 'Jfac ', 'errJ ', 'Jrad']
dfResults = pd.DataFrame(np.concatenate((Rmeso, nuvec, angBN[...,None],\
                np.log10(condnum[...,None]), Jfac.values[...,None], \
                errJ.values[...,None], Jrad.values[...,None]), axis=1), \
                columns=namecol, index=dt4)
with open(fname_out, 'w') as file:
    file.write('# Swarm sats: ' + str(sats) + '\n')
    file.write('# Time-shift in sec.: '+ str(tshift) + '\n')
    file.write('# Along track separation in sec.: '+ str(dt_along) + '\n')    
    file.write('# Use filter: ' + str(use_filter) + '\n')    
    file.write('# time [%Y-%m-%dT%H:%M:%S.%f],\tRc{xyz} [km],\t nuvec{xyz},\t angBN [deg],\
    log_CN,\t  Jfac [microA/m^2],\terrJ [microA/m^2],\t  Jrad [microA/m^2]\n')
dfResults.to_csv(fname_out, mode='a', sep=",", date_format='%Y-%m-%dT%H:%M:%S.%f', \
                 float_format='%15.5f', header=False, na_rep='NaN')


#PLOTS THE RESULTS
#====================================
# creates fig and axes objects
fig_size = (8.27,11.69)
fig = plt.figure(figsize=fig_size, frameon=True)
fig.suptitle('FAC density estimate with dual-sat. method\n', size='xx-large', y=0.995)

# PANEL PART
# designes the panel frames 
xle, xri, ybo, yto = 0.125, 0.95, 0.34, 0.07     # plot margins
ratp = np.array([1, 1, 0.6, 0.6, 1.5, 0.6 ]) # relative hight of each panel 
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
    ax[ii].set_xlim(pd.Timestamp(dtime_beg), pd.Timestamp(dtime_end))

for ii in range(nrp -1):
    ax[ii].set_xticklabels([])
    
#Plots time-series quantities    

ax[0].set_title('Time interval: ' + dtime_beg[:10]+'  ' +\
              dtime_beg[11:19] + ' - '+dtime_end[11:19] +'\n' + \
              'sats: [' + sats[0]+', ' + sats[1]+']        ' + \
              'time-shifts : '+ str(tshift) + ' sec.'\
              '        use_filter: ' + str(use_filter))
ax[0].plot(RsphBswBmod[('Bsw',sats[indmax])] - RsphBswBmod[('Bmod',sats[indmax])])
ax[0].set_ylabel('$dB_{NEC}$ sw'+sats[indmax]+'\n[nT]', linespacing=1.7)
ax[0].legend(['dB_N', 'dB_E', 'dB_C' ], loc = (0.95, 0.1), handlelength=1)

ax[1].plot(RsphBswBmod[('Bsw',sats[indmin])] - RsphBswBmod[('Bmod',sats[indmin])])
ax[1].get_shared_y_axes().join(ax[0], ax[1])
ax[1].get_shared_x_axes().join(ax[0], ax[1])
ax[1].set_ylabel('$dB_{NEC}$ sw'+sats[indmin]+'\n[nT]', linespacing=1.7)
ax[1].legend(['dB_N', 'dB_E', 'dB_C' ], loc = (0.95, 0.1), handlelength=1)

ax[2].plot(dt4, np.log10(condnum))
ax[2].set_ylabel('log(CN)\n')

ax[3].plot(dt4, angBN)
ax[3].set_ylabel('angBN\n[deg]')
if len(ANGind):
    for ii in range(len(grpANGind)):
        ax[3].axvline(ANGstart[ii], ls='--', c='k', lw=0.5)
        ax[3].axvline(ANGstop[ii], ls='--', c='k', lw=0.5)    

ax[4].plot(Jfac)
ax[4].plot(FAC_L2['FAC'])
ax[4].axhline(y=0, linestyle='--', color='k', linewidth=0.7)
ax[4].set_ylabel('$J_{FAC}$\n[$\mu A/m^2$]', linespacing=1.7)
ax[4].legend(['$\mathrm{LS}$', '$\mathrm{ESA\,\, L2}$'], loc = (0.93, 0.1), \
      handlelength=1)

ax[5].plot(errJ)
ax[5].set_ylabel('errJ\n[$\mu A/m^2$]', linespacing=1.7)
if len(EZind):
    for ii in range(len(grpEZind)):
        ax[5].axvline(EZstart[ii], ls='--', c='k', lw=0.5)
        ax[5].axvline(EZstop[ii], ls='--', c='k', lw=0.5)  

# # designes the multiple line xticklabels

locx = ax[5].get_xticks()

dt_trail = dt[:ndt4]
latc = Rsph[:ndt4,indmax,0]
lonc = Rsph[:ndt4,indmax,1]
QDlatc = Aux[:ndt4,indmax,0]
QDlonc = Aux[:ndt4,indmax,1]
MLTc = Aux[:ndt4,indmax,2]

latc_ipl = np.round(np.interp(locx, mdt.date2num(dt_trail), \
                            latc), decimals=2).astype('str')
lonc_ipl = np.round(np.interp(locx, mdt.date2num(dt_trail), \
                            lonc), decimals=2).astype('str')
QDlatc_ipl = np.round(np.interp(locx, mdt.date2num(dt_trail), \
                            QDlatc), decimals=2).astype('str')
QDlonc_ipl = np.round(np.interp(locx, mdt.date2num(dt_trail), \
                            QDlonc), decimals=2).astype('str')
MLTc_ipl = np.round(np.interp(locx, mdt.date2num(dt_trail), \
                            MLTc), decimals=2).astype('str')

lab_fin = ['']*len(locx)
for ii in range(len(locx)):
    lab_ini = mdt.num2date(locx[ii]).strftime('%H:%M:%S')
    lab_fin[ii] = lab_ini + '\n' +latc_ipl[ii] + '\n' + lonc_ipl[ii]+ \
            '\n' + QDlatc_ipl[ii] + '\n' + QDlonc_ipl[ii] + '\n' + MLTc_ipl[ii]
ax[5].set_xticklabels(lab_fin)
plt.figtext(0.01, 0.252, 'Time\nLat\nLon\nQDLat\nQDLon\nMLT')


# INSET PART
# designes the insets to plot s/c constellation 
xmar, ymar, xsep = 0.1, 0.04, 0.05
xwi = (1 - 2*xmar - 2*xsep)/3. # inset width
rat = fig_size[0]/fig_size[1]       # useful to fix the same scales on x and y 

# creates axes for each panel
ax_conf = [0, 0, 0.]
for ii in range(3):
    ax_conf[ii] =  fig.add_axes([xmar + ii*(xwi+xsep), ymar, xwi, xwi*rat])


ic=np.array([1, Rmeso.shape[0]//2, Rmeso.shape[0]-2])   # indexes for ploting the s/c
        
Ri_nec = np.full((len(ic),4,3),np.nan)
Vi_nec = np.full((len(ic),4,3),np.nan)
nlim_nec = np.zeros((len(ic), 2))
elim_nec = np.zeros((len(ic), 2))

z_geo = np.array([0., 0., 1.])
for ii in range(len(ic)):
    pos_i = ic[ii]
    # computes NEC associated with the mesocenter location Rmeso 
    Ci_unit = -Rmeso[pos_i,:]/np.linalg.norm(Rmeso[pos_i,:])[...,None]
    Ei_geo = np.cross(Ci_unit, z_geo)
    Ei_unit = Ei_geo/np.linalg.norm(Ei_geo)
    Ni_geo = np.cross(Ei_unit, Ci_unit)
    Ni_unit = Ni_geo/np.linalg.norm(Ni_geo)
    trmat_i = np.stack((Ni_unit, Ei_unit, Ci_unit))
    # transform the sats. position from GEO mesocentric to NEC mesocentric
    Ri_geo = R4s[pos_i, :, :]
    Ri_nec[ii, :, :] = np.matmul(trmat_i, Ri_geo[...,None]).reshape(Ri_geo.shape)
    nlim_nec[ii,:] = [max(Ri_nec[ii, :, 0]), min(Ri_nec[ii, :, 0])]
    elim_nec[ii,:] = [max(Ri_nec[ii, :, 1]), min(Ri_nec[ii, :, 1])]  
    Vi_geo = (R4s[pos_i+1, :, :] - R4s[pos_i-1, :, :])/2.
    Vi_geo_unit = Vi_geo/np.linalg.norm(Vi_geo,axis=-1)[...,None]
    Vi_nec[ii, :, :] = np.matmul(trmat_i, Vi_geo_unit[...,None]).reshape(Vi_geo.shape)
    
    
# computes the (common) range along N and E
dn_span = max(nlim_nec[:,0] - nlim_nec[:,1])
de_span = max(elim_nec[:,0] - elim_nec[:,1])
d_span = max(np.concatenate((nlim_nec[:,0] - nlim_nec[:,1], 
                             elim_nec[:,0] - elim_nec[:,1])))*1.2

# plots the s/c positions
icolor = ['b', 'r']
for ii in range(len(ic)):
    norig = np.mean(nlim_nec[ii, :])
    eorig = np.mean(elim_nec[ii, :])    
    ax_conf[ii].set_xlim( norig - d_span/2., norig + d_span/2. )
    ax_conf[ii].set_ylim( eorig - d_span/2., eorig + d_span/2. )   
    xquad, yquad = [], []
    for kk in [0, 1, 3, 2, 0]:
        xquad.append(Ri_nec[ii, kk, 0])
        yquad.append(Ri_nec[ii, kk, 1])       
    ax_conf[ii].plot(xquad, yquad, c='k', linestyle=':', linewidth=1)        
    ax_conf[ii].arrow(0, 0, d_span/10*Vi_nec[ii, kk, 0], d_span/10*Vi_nec[ii, kk, 1], \
               color='k', head_width = 4)
    for jj in range(4):
        ax_conf[ii].scatter(Ri_nec[ii, jj, 0], Ri_nec[ii, jj, 1],  marker='o'  ,\
           c=icolor[jj % 2], label=sats[jj%2])      
         
        
ax_conf[0].set_ylabel('East [km]')
ax_conf[0].set_xlabel('North [km]')
ax_conf[1].set_xlabel('North [km]')
ax_conf[1].set_yticklabels([])
ax_conf[2].set_xlabel('North [km]')
ax_conf[2].set_yticklabels([])
handles, labels = ax_conf[2].get_legend_handles_labels()
ax_conf[2].legend(handles[0:2],['sw'+sats[0],'sw'+sats[1]], loc = (0.98, 0.7),\
                  labelcolor=icolor, handlelength=1, ncol = 1)

plt.show()
fig.savefig(fname_fig)


