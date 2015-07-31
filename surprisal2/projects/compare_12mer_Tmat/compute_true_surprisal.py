import numpy as np
import os, sys, glob

import numpy as np
from scipy import loadtxt
from scipy.io import mmread

from msmbuilder.MSMLib import *

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

# default settings
# VAV:  For some reason, this needs to be set in the main body of the script, not in a subroutine
plt.rc('figure', figsize=(3.3, 2.5))
plt.rc('lines', linewidth=1.5, markersize=5)
plt.rc('font', size=12.0)
fontfamily = {'family':'sans-serif','sans-serif':['Arial']}
plt.rc('font', **fontfamily)
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
plt.rc('legend', fontsize=12.0)


sys.path.append('../')
from Figure import Figure

def get_pops(T):
    evals, evecs = GetEigenvectors(T, 6)
    return np.real(evecs[:,0]/evecs[:,0].sum())

def H(row):
    result = 0.
    for j in range(row.shape[0]):
        if row[i] > 0.:
            result += -row[i]*np.log(row[i])
    return result 


if __name__ == '__main__':

  if (1):

    print 'Reading in transition matrices...'
    T_wt = np.array(mmread('../build_12mer_Tmat/macroT_AAAAAA.mtx').todense())
    T_mut = np.array(mmread('../build_12mer_Tmat/macroT_AAFAAA.mtx').todense())
    T_comb = 0.5*(T_wt + T_mut)

    print 'Getting equilibrium populations...'
    pops_wt = get_pops(T_wt)
    pops_mut = get_pops(T_mut)
    pops_comb = get_pops(T_comb)

    print 'T_wt', T_wt
    print 'T_mut', T_mut
    print 'pops_wt', pops_wt
    print 'pops_mut', pops_mut

    surprisals = []
    for i in range(T_wt.shape[0]):
        H_i = H(T_wt[i,:])
        H_i_star = H(T_mut[i,:])
        H_comb = H(T_comb[i,:])
        surprisals.append( H_comb - 0.5*H_i - 0.5*H_i_star )
    surprisals = np.array(surprisals)

    # Harmonic mean
    #pops = pops_wt*pops_mut*2.0/(pops_wt+pops_mut)
    # Average
    pops = (pops_wt+pops_mut)/2.0
    

    # Create the figure
    f = Figure(plt)
    f.add_panel(0.62, 0.15, 0.36, 0.82)     # [xorig, yorig, width, height]
    f.add_panel(0.12, 0.15, 0.36, 0.3)
    f.add_panel(0.12, 0.6, 0.36, 0.37)

   
    ################################
    # Panel A -  s_i versus pi_i
  if (1):
    ax = plt.axes(f.axis_handles[0])

    # plot contours
    sc = np.array([10**i for i in np.arange(-7,-2,0.01)])
    D_values = [10**i for i in range(-8,-4)]
    for D in D_values:
        plt.plot(D/sc, sc, '-', linewidth=1, markersize=3)
    plt.ylim(1e-7, 1e-2)
    plt.yscale('log')
    plt.xlim(1e-4, 1.0)
    plt.xscale('log')
    #plt.legend(D_values)

    # plot pi_i, s_i data
    Ind = ~np.isnan(surprisals)
    plt.plot(pops[Ind], surprisals[Ind], 'ko', markersize=3, markerfacecolor='None')
    UseLabels = True
    if (UseLabels):
      for i in range(len(surprisals)):
        if pops[i]*surprisals[i] > 1e-6:
          print 'State', i, 'pi_i =', pops[i], 's_i =', surprisals[i]
    plt.xlabel('$\\bar{\pi}_i$')
    plt.ylabel('$\\tilde{s}_i$')


    ###############################
    # Panel B - plot pi_i estimate vs. true 
    ax = plt.axes(f.axis_handles[2])
    plt.plot(pops, pops_comb, 'ko', markersize=2, markerfacecolor='None')
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

    CDF = []
    for i in range(len(JS_values)):
       CDF.append( JS_values[0:i+1].sum() )
    print 'JS_values.sum()', JS_values.sum()
    CDF = np.array(CDF)/JS_values.sum()
    print 'CDF', CDF
 
    plt.plot(range(len(JS_values)), CDF, '-')
    plt.xlabel('number of macrostates') 
    plt.ylabel('CDF($\pi_i s_i$)')
    #plt.ylim(1e-7, 1e-2)
    #plt.yscale('log')
    plt.xlim(0, 72)
    #plt.xscale('log')

    ###############################
    # Panel C - plot pi_i estimate vs. true 
    ax = plt.axes(f.axis_handles[2])
    plt.plot(pops, pops_comb, 'ko', markersize=2, markerfacecolor='None')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('$\\bar{\pi}_i$ (estimate)')
    plt.ylabel('$\pi_i$ (true)')


outfile = 'true_surprisal.eps'
plt.savefig(outfile, format='eps')
print 'Wrote:', outfile

