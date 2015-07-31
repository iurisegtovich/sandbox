import numpy as np
import os, sys, glob

import numpy as np
from scipy import loadtxt
from scipy.io import mmread

sys.path.append('../../../src')
from Sampler import *
from Surprisal import *

from msmbuilder.MSMLib import *
from msmbuilder.MSMLib import msm_analysis

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

from Figure import Figure


# default settings
# VAV:  For some reason, this needs to be set in the main body of the script, not in a subroutine
plt.rc('figure', figsize=(2.5, 1.5))
plt.rc('lines', linewidth=1.5, markersize=3)
plt.rc('font', size=12.0)
fontfamily = {'family':'sans-serif','sans-serif':['Arial']}
plt.rc('font', **fontfamily)
plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12)
plt.rc('legend', fontsize=12)




def get_equil_pops(T):
    """Return the stationary eigenvector of a transition matrix T, normalized to populations."""
    evals, evecs = msm_analysis.get_eigenvectors(T, min(T.shape[0],10), 0.001, 50)
    return np.real(evecs[:,0]/evecs[:,0].sum())



def H(row):
    result = 0.
    for j in range(row.shape[0]):
        if row[i] > 0.:
            result += -row[i]*np.log(row[i])
    return result 

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



#########################################

if __name__ == '__main__':

  if (1):

    print 'Reading in transition matrices...'   #tCounts.mtx.ERR		tCounts.mtx.RRR
    C_wt = np.array(mmread('tCounts.MSM.mtx').todense()) 
    C_mut = np.array(mmread('tCounts.MaxCal.mtx').todense()) 
    C_comb = 0.5*(C_wt + C_mut)

    print 'Getting equilibrium populations...'
    # Use the self-consistent method to enforce detailed balance
    T_wt = estimate_transition_matrix(C_wt)  # row-normalized estimate
    T_mut = estimate_transition_matrix(C_mut)  # row-normalized estimate
    T_comb = estimate_transition_matrix(C_comb)  # row-normalized estimate
    T_wt, pops_wt = self_consistent_estimation(T_wt)
    T_mut, pops_mut  = self_consistent_estimation(T_mut)
    T_comb, pops_comb  = self_consistent_estimation(T_comb)

    print 'T_wt', T_wt
    print 'T_mut', T_mut
    print 'pops_wt', pops_wt
    print 'pops_mut', pops_mut

    surprisals = []
    var_surprisals = []
    jsds = []
    var_jsds = []

    nstates = pops_wt.shape[0]
    print 'nstates', nstates

    for k in range(nstates):

        surprisals.append( surprisal(C_wt[k,:], C_mut[k,:], normalize=True) )
        if surprisals[-1] == 0.0:
            print k, 'has s_i=0!  C_wt[k,:], C_mut[k,:] =', C_wt[k,:], C_mut[k,:]
        var_surprisals.append( surprisal_var(C_wt[k,:], C_mut[k,:], normalize=True, bootstrap=False) ) 

        # approximate JSD using pi_si! 8/2014  it's fast but will it work?
        jsds.append( pops_comb[k]*surprisals[-1] )
        var_jsds.append( pops_comb[k]*var_surprisals[-1] )


    # convert quantities to np arrays for ease of use
    surprisals = np.array(surprisals)
    var_surprisals =  np.array(var_surprisals)
    jsds = np.array(jsds)
    var_jsds = np.array(var_jsds)

    # Harmonic mean
    #pops = pops_wt*pops_mut*2.0/(pops_wt+pops_mut)
    # Average
    pops = (pops_wt+pops_mut)/2.0


    # Create the figure
    f = Figure(plt, figsize=(3.25, 2.5), linewidth=1.5, markersize=3,
                       fontsize=9.0, legend_fontsize='medium', xlabelsize = 'small', ylabelsize='small' )
    f.add_panel(0.6, 0.15, 0.38, 0.82)     # [xorig, yorig, width, height]
    f.add_panel(0.12, 0.15, 0.36, 0.3)
    f.add_panel(0.12, 0.6, 0.36, 0.37)

    ################################
    # Panel A -  s_i versus pi_i
  if (1):
    ax = plt.axes(f.axis_handles[0])

    # plot contours
    sc = np.array([10**i for i in np.arange(-8,0,0.01)])
    D_values = [10**i for i in np.arange(-8,0,1)]
    for D in D_values:
        plt.plot(D/sc, sc, '-', linewidth=1)
    plt.xlim(1e-5, 1.0)
    plt.xscale('log')
    plt.ylim(1e-3, 1.0)
    plt.yscale('log')
    #plt.legend(D_values)

    # plot pi_i, s_i data
    Ind = ~np.isnan(surprisals)
    #plt.plot(pops[Ind], surprisals[Ind], 'ko', markerfacecolor='None')
    plt.errorbar(pops[Ind], surprisals[Ind], yerr=var_surprisals[Ind]**0.5, fmt='ko', markerfacecolor='None', elinewidth=0.5, capsize=1)

    UseLabels = True
    if (UseLabels):
      for i in range(len(surprisals)):
        if pops[i]*surprisals[i] > 5e-4:
          print 'State', i, 'pi_i =', pops[i], 's_i =', surprisals[i], 'pi_i*si_i =', pops[i]*surprisals[i]
          plt.text(pops[i], surprisals[i], str(i), fontsize=6, color='r', weight='bold')
    plt.xlabel('$\\bar{\pi}_i$')
    plt.ylabel('$s_i$')


    ###############################
    # Panel B - plot pi_i estimate vs. true 
    ax = plt.axes(f.axis_handles[2])
    plt.plot(pops, pops_comb, 'ko', markerfacecolor='None')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$\\bar{\pi}_i$ (estimate)')
    plt.ylabel('$\pi_i$ (true)')


    ###############################
    # Panel B - Plot the CDF
    ax = plt.axes(f.axis_handles[1])
    JS_values = -pops*surprisals
    JS_values.sort()
    JS_values = -JS_values

    CDF = [0]
    for i in range(len(JS_values)):
       CDF.append( JS_values[0:i+1].sum() )
    print 'JS_values.sum()', JS_values.sum()
    CDF = np.array(CDF)/JS_values.sum()
    print 'CDF', CDF
 
    plt.plot(range(len(JS_values)+1), CDF, '-')
    plt.xlabel('number of macrostates') 
    plt.ylabel('CDF($\pi_i s_i$)')
    #plt.ylim(1e-7, 1e-2)
    #plt.yscale('log')
    plt.xlim(0, nstates)
    #plt.xscale('log')

    ###############################
    # Panel C - plot pi_i estimate vs. true 
    ax = plt.axes(f.axis_handles[2])
    plt.plot(pops, pops_comb, 'ko', markerfacecolor='None')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$\\bar{\pi}_i$ (estimate)')
    plt.ylabel('$\pi_i$ (true)')


outfile = 'fspeptide_RRR_ERR_jsds.eps'

fig = plt.gcf()
fig.set_size_inches(3.5,2.5)
plt.savefig(outfile, format='eps')
print 'Wrote:', outfile

