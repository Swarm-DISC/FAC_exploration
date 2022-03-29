#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 20 12:05:43 2020

@author: blagau
"""

import numpy as np
import matplotlib.patches as mpatches

#SAVES INPUT AND OUTPUT TO ASCII FILES
#====================================
str_trange = dtime_beg.replace('-','')[:8]+'_'+ \
     dtime_beg.replace(':','')[11:17] + '_'+dtime_end.replace(':','')[11:17]

fname_in = 'input_3sat_atmuso_' + str_trange +'.dat'
fname_out = 'FAC_swABC_atmuso_'+ str_trange +'.dat'
fname_FlagsB = 'FlagsB_nonzero_'+ str_trange +'.dat'
fname_fig = 'FAC_swABC_atmuso_'+ str_trange +'.pdf'
# exports the input data
with open(fname_in, 'w') as file:
    file.write('# time [%Y-%m-%dT%H:%M:%S],\tRsph'+str(sats)+' {lat [deg], lon [deg], Radius [m]}, \
    Bsw'+str(sats)+' {NEC [nT]},\tBmod'+str(sats)+' {NEC [nT]}\n')
RsphBswBmod.to_csv(fname_in, mode='a', sep=",", date_format='%Y-%m-%dT%H:%M:%S', \
                   float_format='%15.5f', header=False)
# exports the results
namecol =  ['Rcx','Rcy','Rcz','ux ','uy ','uz ', 'angBN ', 'log_CN3 ', 'Jfac ', 'errJ']
dfResults = pd.DataFrame(np.concatenate((Rc, nuvec, angBN.values[...,None],CN3.values[...,None], 
                     Jfac.values[...,None], errJ.values[...,None]), axis=1),
                   columns=namecol, index=dt).dropna()
with open(fname_out, 'w') as file:
    file.write('# Swarm sats: ' + str(sats) + '\n')
    file.write('# Model: ' + Bmodel + '\n')
    file.write('# Time-shift in sec.: '+ str(ts) + '\n')
    file.write('# Use filter: ' + str(use_filter) + '\n')    
    file.write('# time [%Y-%m-%dT%H:%M:%S.%f],\tRc{xyz} [km],\t nuvec{xyz},\t angBN [deg],\
    log_CN3,\t  Jfac [microA/m^2],\terrJ [microA/m^2]\n')
dfResults.to_csv(fname_out, mode='a', sep=",", date_format='%Y-%m-%dT%H:%M:%S.%f', float_format='%15.5f', 
                     header=False)
# exports times with poor VFM data quality (if applicable)
if (dfFB.sum(axis=0).sum(0) > 0 or datagaps['A'].size) :
    with open(fname_FlagsB, 'w') as file:
        file.write('#T ime stamps with nonzero FlagsB:\n')
        FBnonzero.to_csv(fname_FlagsB,  mode='a', date_format='%Y-%m-%dT%H:%M:%S', sep="\t")
        file.write('\n')
        if datagaps['A'].size :
            file.write('# Missing data for Swarm A: ' + str(datagaps['A'].values)+'\n')    
        if datagaps['B'].size :
            file.write('# Missing data for Swarm B: ' + str(datagaps['B'].values)+'\n')
        if datagaps['C'].size :
            file.write('# Missing data for Swarm C: ' + str(datagaps['C'].values)+'\n')    


#PLOTS THE RESULTS
#====================================
# creates fig and axes objects
fig_size = (8.27,11.69)
fig = plt.figure(figsize=fig_size, frameon=True)
fig.suptitle('FAC density estimate with three-sat. method\n', size='xx-large', y=0.995)

# PANEL PART
# designes the panel frames 
xle, xri, ybo, yto = 0.125, 0.95, 0.3, 0.07     # plot margins
ratp = np.array([1, 1, 1, 0.7, 0.7, 1.5, 0.7 ]) #relative hight of each panel 
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
ax[0].plot(RsphBswBmod[('Bsw','A')] - RsphBswBmod[('Bmod','A')])
ax[0].set_title('Time interval: ' + dtime_beg[:10]+'  ' +\
              dtime_beg[11:19] + ' - '+dtime_end[11:19] +'\n' + \
             'time-shifts [s]: '+ str(ts) + '   use_filter: ' + str(use_filter))
ax[0].set_ylabel('$dB_{NEC}$ swA\n[nT]', linespacing=1.7)
ax[0].legend(['dB_N', 'dB_E', 'dB_C' ], loc = (0.95, 0.1), handlelength=1)

ax[1].plot(RsphBswBmod[('Bsw','B')] - RsphBswBmod[('Bmod','B')])
ax[1].get_shared_y_axes().join(ax[0], ax[1], ax[2])
ax[1].get_shared_x_axes().join(ax[0], ax[1], ax[2])
ax[1].set_ylabel('$dB_{NEC}$ swB\n[nT]', linespacing=1.7)

ax[2].plot(RsphBswBmod[('Bsw','C')] - RsphBswBmod[('Bmod','C')])
ax[2].get_shared_y_axes().join(ax[0], ax[1], ax[2])
ax[2].get_shared_x_axes().join(ax[0], ax[1], ax[2])
ax[2].set_ylabel('$dB_{NEC}$ swC\n[nT]', linespacing=1.7)

ax[3].plot(CN3)
ax[3].set_ylabel('log(CN3)\n')

ax[4].plot(angBN)
ax[4].set_ylabel('angBN\n[deg]\n', linespacing=1.7)

ax[5].plot(Jfac)
ax[5].plot(FAC_L2['FAC'])
ax[5].axhline(y=0, linestyle='--', color='k', linewidth=0.7)
ax[5].set_ylabel('$J_{FAC}$\n[$\mu A/m^2$]', linespacing=1.7)
ax[5].legend(['$\mathrm{J_{ABC}}$', '$\mathrm{J_{L2}}$'], loc = (0.95, 0.1), \
      handlelength=1)

ax[6].plot(errJ)
ax[6].set_ylabel('errJ\n[$\mu A/m^2$]', linespacing=1.7)

# designes the multiple line xticklabels
XYsq = Rc[:,0]**2 + Rc[:,1]**2
Rcr = np.sqrt(XYsq + Rc[:,2]**2)            
latc = np.degrees(np.arctan2(Rc[:,2],np.sqrt(XYsq)))    
lonc = np.degrees(np.arctan2(Rc[:,1], Rc[:,0]))

locx = ax[6].get_xticks()
latc_ipl = np.round(np.interp(locx, mdt.date2num(dt.values), \
                            latc), decimals=2).astype('str')
lonc_ipl = np.round(np.interp(locx, mdt.date2num(dt.values), \
                            lonc), decimals=2).astype('str')
lab_fin = ['']*len(locx)
for ii in range(len(locx)):
    lab_ini = mdt.num2date(locx[ii]).strftime('%H:%M:%S')
    lab_fin[ii] = lab_ini + '\n' +latc_ipl[ii] + '\n' + lonc_ipl[ii]
ax[6].set_xticklabels(lab_fin)
plt.figtext(0.01, 0.255, 'Time\nLat\nLon')


# INSET PART
# designes the insets to plot s/c constellation 
xmar, ymar, xsep = 0.1, 0.04, 0.05
xwi = (1 - 2*xmar - 2*xsep)/3. # inset width
rat = fig_size[0]/fig_size[1]       # useful to fix the same scales on x and y 

# creates axes for each panel
ax_conf = [0, 0, 0.]
for ii in range(3):
    ax_conf[ii] =  fig.add_axes([xmar + ii*(xwi+xsep), ymar, xwi, xwi*rat])


ic=np.array([1, R.shape[0]//2, R.shape[0]-2])   # indexes for ploting the s/c
# for ii in range(len(ax)):                         # uncomment to plot vertical lines at dt[ic]
#     for jj in range(len(ic)):
#         ax[ii].axvline(x=dt[ic[jj]], ls=(0,(5,10)), lw=1, color='k')
        
Ri_nec = np.full((len(ic),nsc,3),np.nan)
Vi_nec = np.full((len(ic),nsc,3),np.nan)
nlim_nec = np.zeros((len(ic), 2))
elim_nec = np.zeros((len(ic), 2))

z_geo = np.array([0., 0., 1.])
for ii in range(len(ic)):
    pos_i = ic[ii]
    # computes NEC associated with the mesocenter location Rc 
    Ci_unit = -Rc[pos_i,:]/np.linalg.norm(Rc[pos_i,:])[...,None]
    Ei_geo = np.cross(Ci_unit, z_geo)
    Ei_unit = Ei_geo/np.linalg.norm(Ei_geo)
    Ni_geo = np.cross(Ei_unit, Ci_unit)
    Ni_unit = Ni_geo/np.linalg.norm(Ni_geo)
    trmat_i = np.stack((Ni_unit, Ei_unit, Ci_unit))
    # transform the sats. position from GEO mesocentric to NEC mesocentric
    Ri_geo = Rmeso[pos_i, :, :]
    Ri_nec[ii, :, :] = np.matmul(trmat_i, Ri_geo[...,None]).reshape(Ri_geo.shape)
    nlim_nec[ii,:] = [max(Ri_nec[ii, :, 0]), min(Ri_nec[ii, :, 0])]
    elim_nec[ii,:] = [max(Ri_nec[ii, :, 1]), min(Ri_nec[ii, :, 1])]  
    Vi_geo = (R[pos_i+1, :, :] - R[pos_i-1, :, :])/2.
    Vi_geo_unit = Vi_geo/np.linalg.norm(Vi_geo,axis=-1)
    Vi_nec[ii, :, :] = np.matmul(trmat_i, Vi_geo_unit[...,None]).reshape(Vi_geo.shape)
    
    
# computes the (common) range along N and E
dn_span = max(nlim_nec[:,0] - nlim_nec[:,1])
de_span = max(elim_nec[:,0] - elim_nec[:,1])
d_span = max(np.concatenate((nlim_nec[:,0] - nlim_nec[:,1], 
                             elim_nec[:,0] - elim_nec[:,1])))*1.2

# plots the s/c positions
icolor = ['b', 'g', 'r']
for ii in range(len(ic)):
    norig = np.mean(nlim_nec[ii, :])
    eorig = np.mean(elim_nec[ii, :])    
    ax_conf[ii].set_xlim( norig - d_span/2., norig + d_span/2. )
    ax_conf[ii].set_ylim( eorig - d_span/2., eorig + d_span/2. )   
    for jj in range(3):
        ax_conf[ii].scatter(Ri_nec[ii, jj, 0], Ri_nec[ii, jj, 1], \
           c=icolor[jj], label = icolor[jj])
        ax_conf[ii].arrow(Ri_nec[ii, jj, 0] , Ri_nec[ii, jj, 1], \
               d_span/10*Vi_nec[ii, jj, 0], d_span/10*Vi_nec[ii, jj, 1], \
               color=icolor[jj], head_width = 4)

ax_conf[0].set_ylabel('East [km]')
ax_conf[0].set_xlabel('North [km]')
ax_conf[1].set_xlabel('North [km]')
ax_conf[1].set_yticklabels([])
ax_conf[2].set_xlabel('North [km]')
ax_conf[2].set_yticklabels([])
label_patches = [
    mpatches.Patch(color=c, label=l) for c, l in zip(icolor, ('swA', 'swB', 'swC'))
]
ax_conf[2].legend(
    handles=label_patches, loc=(0.98, 0.7), handlelength=1, ncol = 1
)


plt.show()

fig.savefig(fname_fig)
