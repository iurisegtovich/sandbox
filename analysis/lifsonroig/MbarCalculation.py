import numpy as np
import pymbar
import mbarTools_GZ as mbartools


ntemps = 24
mintemp = 300.
maxtemp = 800.
dlogtemp = (np.log(maxtemp) - np.log(mintemp))/(ntemps-1)
temps = np.exp(np.arange(np.log(mintemp), np.log(maxtemp)+dlogtemp, dlogtemp ))
temps = np.round(temps)[0:ntemps]

K = ntemps
N = 1000

energies_LogL = np.zeros((K,N),np.float64)
obs_w = np.zeros((K,N),np.float64)
obs_v = np.zeros((K,N),np.float64)
bins=1
for seqlength in range(3,16):
    for k in range(ntemps):
        path_dir = "SPE%d_trj%d_cis-_gaff_u"%(seqlength,k)
        path_LogL = os.path.join(path_dir,'LogLikelihood.dat')
        path_w = os.path.join(path_dir,'w_params.dat')
        path_v = os.path.join(path_dir,'v_params.dat')
        energies[k,:] = np.loadtxt(path_LogL)[-N:]
        obs_w[k,:] = np.loadtxt(path_w)[-N:]
        obs_v[k,:] = np.loadtxt(path_v)[-N:]

mbar = mbartools.ComputeFreeEnergiesFromREMD(energies=energies_LogL, obs=obs_w, bins=1,
                                             temps=temps, N_max = 1000, NSkipTraj=1)


