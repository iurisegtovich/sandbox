import os,sys
import numpy as np
from pymbar import MBAR
import mbarTools_GZ as mbartools


def prepare_mbar_input(u_k, temps, N_max = 1000):

    kB_in_kJ = 1.381e-23 * 6.022e23 / 1000.0 # Boltzmann constant in kJ/mol/K -- for GROMACS
    kB = kB_in_kJ/4.1868   # use for AMBER (kcal/mol/K)
    print 'kB*T =', kB*300., 'kcal/mol at 300K'   # testing

    ntemps = len(temps) # number of temperatures
    K = ntemps # ntemps

    energies = energies[:,::NSkipTraj]

    # Allocate storage for simulation data
    N_k = np.zeros([K], np.int32) # N_k[k] is the number of snapshots from umbrella simulation k
    beta_k = np.zeros([K], np.float64) # beta_k[k] is the inverse temperature for replica k
    u_kn = np.zeros([K,N_max], np.float64) # u_kn[k,n] is the reduced potential energy without umbrella restraints of snapshot n of umbrella simulation k
    u_kln = np.zeros([K,K,N_max], np.float64) # u_kln[k,l,n] is the reduced potential energy of snapshot n from umbrella simulation k evaluated at umbrella l

    # Set the temperatures
    for k in range(ntemps):
        # Parse line k.
        beta_k[k] = 1.0 / (kB * temps[k])
    print 'beta_k', beta_k

    # Read in the energy data
    for k in range(ntemps):
        N_k[k] = N_max
        for n in range(N_k[k]):
            u_kn[k,n] = energies[k,n]
            for l in range(ntemps):
                u_kln[k,l,n] = beta_k[l] / beta_k[k] * u_kn[k,n]
    return u_kln, N_k


ntemps = 24
mintemp = 300.
maxtemp = 800.
dlogtemp = (np.log(maxtemp) - np.log(mintemp))/(ntemps-1)
temps = np.exp(np.arange(np.log(mintemp), np.log(maxtemp)+dlogtemp, dlogtemp ))
temps = np.round(temps)[0:ntemps]

K = ntemps
N = 1000

path_data = "/Volumes/Guangfeng/results/Lifsonroig/peptoid/SPE_lifsonroig_cisminus_new/results"

energies_LogL = np.zeros((K,N),np.float64)
obs_w = np.zeros((K,N),np.float64)
obs_v = np.zeros((K,N),np.float64)
bins=1
for seqlength in range(5,6):
    for k in range(ntemps):
        path_dir = os.path.join(path_data,"SPE%d_trj%d_cis-_gaff_u"%(seqlength,k))
        path_LogL = os.path.join(path_dir,'LogLikelihood.dat')
        path_w = os.path.join(path_dir,'w_params.dat')
        path_v = os.path.join(path_dir,'v_params.dat')
        energies_LogL[k,:] = np.loadtxt(path_LogL)[-N:]
        obs_w[k,:] = np.loadtxt(path_w)[-N:,0]
        obs_v[k,:] = np.loadtxt(path_v)[-N:,0]


u_kln, N_k = prepare_mbar_input(u_k=energies_LogL,temps=temps,N_max=N)

print "Running MBAR..."
mbar = MBAR(u_kln, N_k, relative_tolerance = 2.0e-2, verbose = True, method = 'self-consistent-iteration')
# compute observable expectations
P, dP = mbar.computeExpectations(obs_w,uncertainty_method='approximate')
print "P:"
print P
print "dP"
print dP



