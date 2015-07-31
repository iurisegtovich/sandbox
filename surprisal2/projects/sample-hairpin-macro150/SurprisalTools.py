import os, sys, glob, string
import numpy as np

from math import ceil, fabs
from scipy.io import mmread, mmwrite
from scipy import sparse

from msmbuilder.msm_analysis import get_eigenvectors

sys.path.append('../../src')
from Sampler import *
from Surprisal import *

sys.path.append('../../scripts')
from HelperTools import *
from ssaCalculatorVAV import *
from ssaTools import *

from msmbuilder import MSMLib
from msmbuilder import msm_analysis
#from Reversible_MLE_with_Guess import *

UsePyplot = False #True

if UsePyplot:
    import matplotlib
    from matplotlib import pyplot as plt


def self_consistent_estimation(T, pops=None, tol=1e-6, Verbose=False):
    """Given initial estimates of T and pops, iterate to get revised estimates of
       both, that satisfy detailed-balance."""

    if pops == None:
        pops = get_equil_pops(T)

    epsilon = 1. 
    iteration = 0
    while epsilon > tol: 
        new_T = np.identity(T.shape[0], dtype=float)
        for i in range(T.shape[0]):
          for j in range(T.shape[1]):
            if not (i==j):
                # skip if both ij and ji and zero  
                if not (T[i,j] == 0):
                    a_ij = min(1., pops[j]*T[j,i]/(pops[i]*T[i,j]))
                    new_T[i,j] = a_ij*T[i,j]
                    new_T[i,i] -= a_ij*T[i,j]
        new_pops = get_equil_pops(new_T)
        pops = get_equil_pops(T)

        # Calculate the norm of new vs old 
        epsilon = norm(T,new_T)

        T = new_T
        pops = new_pops
        iteration += 1
        if Verbose: print iteration, epsilon, pops, T

    #print 'self-consistent iterations needed =', iteration ## its always 2!  i.e. this is not an intertive algorithm -- just a once-pass
    return T, pops
        
def norm(T1,T2):
    """Return the norm of two transition matrices; average 1/N \sum_i=1^N \sum_j |T1_ij - T2_ij|."""
    return (np.abs(T1 - T2)).sum()/float(T1.shape[0]) 
    

def get_equil_pops(T):
    """Return the stationary eigenvector of a transition matrix T, normalized to populations."""
    evals, evecs = msm_analysis.get_eigenvectors(T, min(T.shape[0],10), 0.001, 50)
    return np.real(evecs[:,0]/evecs[:,0].sum())



def get_populations(count_mtx, method='mle', initial_guess=None, Verbose=False):
    """Return the estimated equilibirum populations, and symmetrized count matrix"""

    if method == 'transpose':
        #counts = count_mtx.todense()
        counts_sym = 0.5*(count_mtx.transpose()+count_mtx)
        #return np.array(counts_sym.sum(axis=1)/counts_sym.sum()).flatten()
        T = MSMLib.estimate_transition_matrix(counts_sym)

    if method == 'mle':
        #print "count_mtx.sum() before MLE", count_mtx.sum()
        try:
            counts_sym = MSMLib.mle_reversible_count_matrix(count_mtx, prior=1.0, initial_guess=initial_guess)
            #counts_sym = MSMLib.mle_reversible_count_matrix(count_mtx, prior=0.0, initial_guess=initial_guess)
        except:
            if Verbose: print "initial guess not working. Using guess-less mle"
            counts_sym = MSMLib.mle_reversible_count_matrix(sparse.csr_matrix(count_mtx))
        #print "counts_sym.sum() after MLE", counts_sym.sum()
        T = estimate_transition_matrix(counts_sym)

    evals, evecs = msm_analysis.get_eigenvectors(T, min(count_mtx.shape[0],10), 0.001, 50)
    return T, np.real(evecs[:,0]/evecs[:,0].sum())  # T and pops


def get_true_populations(tProb):
    evals, evecs = msm_analysis.get_eigenvectors(tProb, min(tProb.shape[0],10), 0.001, 50)
    return np.real(evecs[:,0]/evecs[:,0].sum())


def get_equilibrium_counts(tProb, ncounts):
    """Return a matrix of ncounts total counts, portioned according to equilibrium populations."""

    tProb = np.array( tProb.todense() )

    evals, evecs = msm_analysis.get_eigenvectors(tProb, min(tProb.shape[0],10), 0.001, 50)
    populations =  np.real(evecs[:,0]/evecs[:,0].sum()) 

    counts = np.zeros(tProb.shape) 
    for i in range(len(populations)):
        counts[i,:] = ncounts*populations[i]*tProb[i,:]    
    return counts


def dH_dn(counts):
    """Calculate the derivative of the entropy with respect to the counts n ."""

    counts = counts.astype(float)
    total = counts.sum()
    result = np.zeros( counts.shape[0] ) 
    for i in range(counts.shape[0]):
        result[i] = (total - counts[i])/(total**2) * np.log( (total - counts[i])/counts[i] )
    return result


def estimate_JSD(rowcounts1, rowcounts2, pop1, pop2, pop_comb, UsePseudocounts=True, nbootstraps = 100):
    """Estimate the mean and variance of the JSD_i for row i, given the state populations.
    """

    ncounts1 = float(rowcounts1.sum())  # the total number of counts given to us
    ncounts2 = float(rowcounts2.sum()) 

    if UsePseudocounts:
        padded_counts1 = rowcounts1 + 1.0  #/len(rowcounts1)
        padded_counts2 = rowcounts2 + 1.0  #/len(rowcounts2)
        posterior_T1 = padded_counts1/float(padded_counts1.sum())
        posterior_T2 = padded_counts2/float(padded_counts2.sum())
    else:
        posterior_T1 = rowcounts1/float(rowcounts1.sum())
        posterior_T2 = rowcounts2/float(rowcounts2.sum())

    # Draw a bunch of resampled counts  (if the counts are fractional draw at least *one* count)
    if ncounts1 > 0:
        boot1 = np.random.multinomial(max(1,ncounts1), posterior_T1, size=nbootstraps).astype(float) # returns size (nbootstraps, nstates)
    else:
        boot1 = np.zeros( (nbootstraps, len(rowcounts1)) )
    if ncounts2 > 0:
        boot2 = np.random.multinomial(max(1,ncounts2), posterior_T2, size=nbootstraps).astype(float)
    else:
        boot2 = np.zeros( (nbootstraps, len(rowcounts2)) )

    boot_comb = boot1+boot2
   
    p = ncounts1/(ncounts1+ncounts2)
    naive_pop_comb = p*pop1+(1.-p)*pop2
    JSD_boot = np.array( [ (naive_pop_comb*H(boot_comb[k,:]) - p*pop1*H(boot1[k,:]) - (1.-p)*pop2*H(boot2[k,:])) for k in range(nbootstraps) ] ) 
    return JSD_boot.mean(), JSD_boot.var()


def get_JSD(rowcounts1, rowcounts2, pop1, pop2, pop_comb):
    """Return a non-bootstrapped version of the JSD."""

    counts_comb = rowcounts1 + rowcounts2
    p = rowcounts1.sum()/(rowcounts1.sum()+rowcounts2.sum())
    #return  (pop_comb*H(counts_comb) - pop1*p*H(rowcounts1) - pop2*(1.-p)*H(rowcounts2))
    naive_pop_comb = p*pop1 + (1.-p)*pop2
    return  (naive_pop_comb*H(counts_comb) - pop1*p*H(rowcounts1) - pop2*(1.-p)*H(rowcounts2))



def estimate_surprisal(rowcounts1, rowcounts2, UsePseudocounts=True, nbootstraps = 100):
    """Estimate the mean and variance of the JSD_i for row i, given the state populations.
    """

    ncounts1 = float(rowcounts1.sum())  # the total number of counts given to us
    ncounts2 = float(rowcounts2.sum())

    if UsePseudocounts:
        padded_counts1 = rowcounts1 + 1.0
        padded_counts2 = rowcounts2 + 1.0
        posterior_T1 = padded_counts1/float(padded_counts1.sum())
        posterior_T2 = padded_counts2/float(padded_counts2.sum())
    else:
        posterior_T1 = rowcounts1/float(rowcounts1.sum())
        posterior_T2 = rowcounts2/float(rowcounts2.sum())

    # Draw a bunch of resampled counts - note we have to resample at least one!
    boot1 = np.random.multinomial(max(1,ncounts1), posterior_T1, size=nbootstraps) # returns size (nbootstraps, nstates)
    boot2 = np.random.multinomial(max(1,ncounts2), posterior_T2, size=nbootstraps)
    boot_comb = boot1+boot2
  
    p = ncounts1/(ncounts1+ncounts2)
    JSD_boot = np.array( [ (pop_comb*H(boot_comb[i,:]) - pop1*p*H(boot1[i,:]) - pop2*(1.-p)*H(boot2[i,:])) for i in range(nbootstraps) ] )
    #print 'JSD_boot', JSD_boot

    # JSD_boot = np.array( [ (pop_comb*H(boot_comb[i,:]) - 0.5*pop1*H(boot1[i,:]) - 0.5*pop2*H(boot2[i,:])) for i in range(nbootstraps) ] ) 
    return JSD_boot.mean(), JSD_boot.var()

