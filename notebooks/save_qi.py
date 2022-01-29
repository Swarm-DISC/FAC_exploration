tq_beg = min([qorbs_dB[0][jj].index[0], qorbs_dB[1][jj].index[0]])
tq_end = max([qorbs_dB[0][jj].index[-1], qorbs_dB[1][jj].index[-1]])

str_trange = tq_beg.isoformat().replace('-','')[:8]+'_'+ \
        tq_beg.isoformat().replace(':','')[11:17] + \
        '_'+tq_end.isoformat().replace(':','')[11:17]

with open(fname_out, "a") as file:
    file.write('\n')
    file.write('========\n')
    file.write('Time interval: ' + tq_beg.isoformat()[:10]+'  ' + \
              tq_beg.isoformat()[11:19] + ' - '+tq_end.isoformat()[11:19] +'\n\n')

for ss in range(2):
    str_tmva = tbeg_mva[ss][jj].strftime("%Y-%m-%d   %H:%M:%S") +\
           ' - ' + tend_mva[ss][jj].strftime("%H:%M:%S")
    B_unit = np.squeeze(dir_bunit[ss, jj, :])
    mindir = np.squeeze(dir_min[ss, jj, :])
    maxdir = np.squeeze(dir_max[ss, jj, :]) 
    eigval=[0, lbd_min[ss,jj], lbd_max[ss,jj]]
    
    with open(fname_out, "a") as file:
        file.write('Satellite:  sw' + sats[ss] + '\n')
        file.write('   Orbit nr.: ' + str(int(orbs[ss][0]) + int(jj/4)) + \
              ',  quadrant ' + str(jj - 4*int(jj/4) + 1) +'\n')
        file.write('   MVA interval: ' + str_tmva + '\n')    
        file.write('      B_unit:  '+ format(eigval[0],'.2f')+'  ['+ \
                   format(B_unit[0],'.3f') + ', '+format(B_unit[1], '.3f') +\
                   ', '+format(B_unit[2], '.3f') +'] \n')    
        file.write('      minvar:  '+ format(eigval[1],'.2f')+'  ['+ \
                   format(mindir[0],'.3f') + ', '+format(mindir[1], '.3f') +\
                   ', '+format(mindir[2], '.3f') +'] \n')
        file.write('      maxvar:  '+ format(eigval[2],'.2f')+'  ['+ \
                       format(maxdir[0],'.3f') + ', '+format(maxdir[1], '.3f') +\
                   ', '+format(maxdir[2], '.3f') +'] \n') 
        file.write('   Eigenvalues ratio: ' + format(eigval[2]/eigval[1], '.1f') +'\n')    
        file.write('   FAC inclination wrt sat. velocity (tangential plane):  ' + \
        format(ang_vn[ss,jj], '.1f') + '  deg. \n\n')

with open(fname_out, "a") as file:
    file.write('Correlation analysis: \n')
    file.write('   reference interval: MVA interval on sw'+sats[iref_arr[jj]]+'\n')
    file.write('   correlation coeff.: '+ str(np.round(cc_ls[jj], decimals=3)) + '\n')
    file.write('   optimum time lag [s]: ' + str(opt_lag_ls[jj]) + '\n\n')         
