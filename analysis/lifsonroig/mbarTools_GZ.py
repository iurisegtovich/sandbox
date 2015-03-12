#!/usr/bin/env python

### ADAPTED FROM THE MBAR EXAMPLE SCRIPT: ###
# Example illustrating the application of MBAR to compute a 1D PMF from an umbrella sampling simulation.
#
# The data represents an umbrella sampling simulation for the chi torsion of a valine sidechain in lysozyme L99A with benzene bound in the cavity.
# 
# REFERENCE
# 
# D. L. Mobley, A. P. Graves, J. D. Chodera, A. C. McReynolds, B. K. Shoichet and K. A. Dill, "Predicting absolute ligand binding free energies to a simple model site," Journal of Molecular Biology 371(4):1118-1134 (2007).
# http://dx.doi.org/10.1016/j.jmb.2007.06.002

import os, sys, glob, cPickle

import numpy as np
from math import *
from scipy import loadtxt, savetxt

from pymbar import MBAR # multistate Bennett acceptance ratio
# import pmftools  # tools for making PMFs using MBAR


# Constants.
kB_in_kJ = 1.381e-23 * 6.022e23 / 1000.0 # Boltzmann constant in kJ/mol/K -- for GROMACS
kB = kB_in_kJ/4.1868   # use for AMBER (kcal/mol/K)

print 'kB*T =', kB*300., 'kcal/mol at 300K'   # testing


### Functions ###


def ComputeFreeEnergiesFromREMD(energies, obs, bins, temps, N_max = 500000, NSkipTraj=1): 
    """Computes MBAR estimates F(bins), from

    energies -     a KxN matrix of N energies for trajectories at temperature index K
    obs      -     a KxN matrix of N observable values for trajectories at temperature index K
    temps    -     a length K array of temperatures (in Kelvin) 
    bins     -     a range of left-edge values for binning observable values into states

    Options

    N_max    -     maximum number of snapshots for array allocation (default: 500000)
    NSkipTraj -    number of values to skip when reading in energies/observables (default: 1, i.e. no skip)

    Returns

    P, dP    -     a (nbins)x(K) matrix of equilibrium probability estimates and uncertainties

    """


    # Constants.
    kB_in_kJ = 1.381e-23 * 6.022e23 / 1000.0 # Boltzmann constant in kJ/mol/K -- for GROMACS
    kB = kB_in_kJ/4.1868   # use for AMBER (kcal/mol/K)
    print 'kB*T =', kB*300., 'kcal/mol at 300K'   # testing

    ntemps = len(temps) # number of temperatures
    K = ntemps # ntemps 
    beta = 1.0 / (kB * 300.0) # inverse temperature of simulations (in 1/(kcal/mol))
    #nbins = len(bins)  # number of bins
    #binwidth = bins[1] - bins[0]

    energies = energies[:,::NSkipTraj]
    obs = obs[:,::NSkipTraj]

    nsnaps = energies.shape[1]

    SubsampleData = False

    N_max = min(nsnaps, N_max)

    # Allocate storage for simulation data
    N_k = np.zeros([K], np.int32) # N_k[k] is the number of snapshots from umbrella simulation k
    beta_k = np.zeros([K], np.float64) # beta_k[k] is the inverse temperature for replica k
    bin_kn = np.zeros([K,N_max], np.int32) # bin_kn[k,n] is the PY1-PY@ probe distance for snapshot n from replica simulation k
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
            #u_kn[k,n] = beta_k[k] * energies[k,n]
            for l in range(ntemps):
                u_kln[k,l,n]   = beta_k[l] / beta_k[k] * u_kn[k,n]
    """
            # Compute bin assignment.
            bin_indicators = np.where((obs[k,n]>=bins)&(obs[k,n]<(bins+binwidth)),1,0)
            if np.sum(bin_indicators) > 0:
                bin_kn[k,n] = np.argmax(bin_indicators)
            else:
                bin_kn[k,n] = -1

        print 'bin_kn[k,:]', bin_kn[k,:]

    if SubsampleData:
        #Compute correlation times for potential energy and chi timeseries.
        indices = timeseries.subsampleCorrelatedData(u_kn[k,:])

        # Subsample data.
        N_k[k] = len(indices)
        u_kn[k,0:N_k[k]] = u_kn[k,indices]
        bin_kn[k,0:N_k[k]] = bin_kn[k,indices]

    # Set zero of u_kn -- this is arbitrary.
    u_kn -= u_kn.min()


    print 
    """
    # Initialize MBAR.
    print "Running MBAR..."
    #mbar = MBAR.MBAR(u_kln, N_k, verbose = True, method = 'Newton-Raphson') # doesn't work too well...
    #mbar = MBAR(u_kln, N_k, relative_tolerance = 1.0e-7, verbose = True, method = 'self-consistent-iteration')
    mbar = MBAR(u_kln, N_k, relative_tolerance = 2.0e-2, verbose = True, method = 'self-consistent-iteration')
    return mbar
    """
    print 'binindex\t-binvalue\tP\tdP'
    P, dP = [],[]
    for i in range(nbins):
        A_kn = np.where(bin_kn==i,1,0)
        (p_i, dp_i) = mbar.computeExpectations(A_kn, uncertainty_method='approximate')
        P.append(p_i)
        dP.append(dp_i)
        print i, bins[i], 
        for p in p_i: print p,
        for dp in dp_i: print dp,
        print

    return np.array(P), np.array(dP)
    """
