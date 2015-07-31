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

UsePyplot = True

if UsePyplot:
    import matplotlib
    from matplotlib import pyplot as plt
    matplotlib.use('Agg') 

from SurprisalTools import *


#====================#
# Main

usage = """Usage: python prod_sampling3.py expfile outname

    Example:  python prod_sampling3.py pi_si_variance_uniform_wtmut.exp out"""


if len(sys.argv) < 3:
    print usage
    sys.exit(1)

expfile = sys.argv[1]
outname = sys.argv[2]

from Experiment import *
e = Experiment(expfile)



if (e.parms['system'].count('toy')):
    #====================#
    # toymer hairpin MSM #
    #====================#
  
    # Read in the Transition matrix from file
    wt_tprob = np.array(mmread('../build_12mer_Tmat/macroT_AAAAAA.mtx').todense())
    #mut_tprob = np.array(mmread('../build_12mer_Tmat/macroT_AAAAAA.mtx').todense())  # convergence sampling
    mut_tprob = np.array(mmread('../build_12mer_Tmat/macroT_AAFAAA.mtx').todense())

    wt_initial_counts = np.zeros( wt_tprob.shape )
    mut_initial_counts = np.zeros( mut_tprob.shape )

    # Instantiate a Sampler object
    # -- if the *.jumps file does not exist, it will be calculated and saved to file (as a picked dict)
    wt_model  = Sampler(wt_tprob, 'macroT_AAAAAA.jumps')
    mut_model = Sampler(wt_tprob, 'macroT_AAAAAA.jumps')     # convergence sampling 
    #mut_model = Sampler(mut_tprob, 'macroT_AAFAAA.jumps')   # variance sampling

    wt_true_pops = get_true_populations(wt_tprob)
    mut_true_pops = get_true_populations(mut_tprob)
    comb_true_pops = get_true_populations(wt_tprob + mut_tprob)

    nstates = wt_tprob.shape[0]

    # To start with, sample counts from each state, with a uniform ergodic prior to ensure connectivity
    prior = 0.01
    C_wt  = np.zeros( wt_model.tProb.shape ) + prior*np.ceil(wt_model.tProb)
    C_mut = np.zeros( mut_model.tProb.shape ) + prior*np.ceil(mut_model.tProb)

    n_presamples = 10000

    if (e.parms['system'] == 'toy wt'):
        # MC sample to get initial C_wt
        print 'Sampling %d MC steps for wt...'%n_presamples
        j_wt = np.random.randint(nstates)
        C_wt, last_j_wt  = wt_model.sample_traj(j_wt, n_presamples, previous_counts=C_wt)

    elif (e.parms['system'] == 'toy wt mut'):

        #if e.parms['sampling_method'] == 'uniform':
        # MC sample to get initial C_wt
        print 'Pre-sampling %d MC steps, proportional to population for each state of wt and mut ...'%n_presamples
        j_wt = np.random.randint(nstates)
        j_mut = np.random.randint(nstates)
        for k in range(nstates):
            C_wt  = wt_model.sample(k, int(n_presamples*wt_true_pops[k])+1, previous_counts=C_wt)
            C_mut = wt_model.sample(k, int(n_presamples*mut_true_pops[k])+1, previous_counts=C_mut)
        #else:
        # MC sample to get initial C_wt
        #print 'Pre-sampling %d MC steps for wt...'%n_presamples
        #j_wt = np.random.randint(nstates)
        #C_wt, last_j_wt  = wt_model.sample_traj(j_wt, n_presamples, previous_counts=C_wt)
        #print 'Pre-sampling %d MC steps for mut...'%n_presamples
        #j_mut = np.random.randint(nstates)
        #C_mut, last_j_mut  = mut_model.sample_traj(j_mut, n_presamples, previous_counts=C_mut)



if (e.parms['system'] == 'hairpin-macro150'):
    #=================================#
    # GB1 vs trpzip4 hairpin macro150 #
    #=================================#

    # Read in the Transition matrix from file
    wt_tprob = np.array(mmread('macro150-4x/wt-tProb.mtx').todense())
    mut_tprob = np.array(mmread('macro150-4x/tz4-tProb.mtx').todense())

    # read in initial counts from file
    wt_initial_counts = np.array(mmread('initial-counts/wt-tCounts.mtx').todense())
    mut_initial_counts = np.array(mmread('initial-counts/tz4-tCounts.mtx').todense())

    # Instantiate a Sampler object
    # -- if the *.jumps file does not exist, it will be calculated and saved to file (as a picked dict)
    wt_model  = Sampler(wt_tprob, 'macroT_wt.jumps')
    mut_model = Sampler(wt_tprob, 'macroT_tz4.jumps')

    wt_true_pops = get_true_populations(wt_tprob)
    mut_true_pops = get_true_populations(mut_tprob)
    comb_true_pops = get_true_populations(wt_tprob + mut_tprob)

    nstates = wt_tprob.shape[0]

    # Use initial counts,  with a uniform ergodic prior to ensure connectivity
    scaling = 0.001
    prior = 0.   #0.01 
    C_wt  = scaling*wt_initial_counts + prior*np.ceil(wt_model.tProb)
    C_mut = scaling*mut_initial_counts + prior*np.ceil(mut_model.tProb)

    n_presamples = e.parms['npresamples']

    if n_presamples > 0:
        # MC sample to get initial C_wt
        print 'Pre-sampling %d MC steps, proportional to population for each state of wt and mut ...'%n_presamples
        j_wt = np.random.randint(nstates)
        j_mut = np.random.randint(nstates)
        for k in range(nstates):
            C_wt  = wt_model.sample(k, int(n_presamples*wt_true_pops[k])+1, previous_counts=C_wt)
            C_mut = wt_model.sample(k, int(n_presamples*mut_true_pops[k])+1, previous_counts=C_mut)


    # pick the native state as the first starting state 
    #j = np.random.randint(nstates)
    j = 66 # the native state
    j_wt = j
    j_mut = j


#==========================
# output files
fout_surprisals = open(outname+'.surprisals.dat', 'w')
fout_jsds = open(outname+'.jsds.dat', 'w')
fout_pops = open(outname+'.pops.dat', 'w')

# make some nice headers for printing output

jsd_headers  = ['stateWT', 'stateMUT' 'nsamples', 'JSD', 'var(JSD)']
for k in range(nstates):
    jsd_headers.append( 'JSD_'+str(k) )
for k in range(nstates):
    jsd_headers.append( 'var(JSD)_'+str(k) )
fout_jsds.write( '#'+string.joinfields(jsd_headers, '\t')+'\n' )

surprisal_headers  = ['stateWT', 'stateMUT', 'nsamples', 's', 'var(s)']
for k in range(nstates):
    surprisal_headers.append( 's_'+str(k) )
for k in range(nstates):
    surprisal_headers.append( 'var(s)_'+str(k) )
fout_surprisals.write( '#'+string.joinfields(surprisal_headers, '\t')+'\n' )

pops_headers  = ['stateWT', 'stateMUT', 'nsamples']
for k in range(nstates):
    pops_headers.append( 'pop_'+str(k) )
fout_pops.write( '#'+string.joinfields(pops_headers, '\t')+'\n' )


# sample trajectories
nrounds = e.parms['nrounds']
nsamples = e.parms['nsamples']
total_samples = 0
total_JSD, total_JSD_var = [], []

#===========================================================================#
#sampling_method = 'traj' # 'uniform', 'traj'
#adaptive_method = 'continue traj'  #'JSD', 'JSD variance', 'random', 'continue traj', 'surprisal', 'surprisal variance', 'pi_si', 'pi_si variance'
sampling_method = e.parms['sampling_method']
adaptive_method = e.parms['adaptive_method']
print 'sampling_method', sampling_method
print 'adaptive_method', adaptive_method
#===========================================================================#

for i in range(nrounds):

    if sampling_method == 'uniform':
        C_wt   = wt_model.sample(j_wt, nsamples, previous_counts=C_wt)
        C_mut  = mut_model.sample(j_mut, nsamples, previous_counts=C_mut)

    if sampling_method == 'traj':
        C_wt, final_state_wt   = wt_model.sample_traj(j_wt, nsamples, previous_counts=C_wt)
        C_mut, final_state_mut = mut_model.sample_traj(j_mut, nsamples, previous_counts=C_mut)

    C_comb = C_wt + C_mut

    total_samples += nsamples
    print 'round', i, 'total samples =', total_samples

    # First, we estimate the equilibirum populations
    UseSelfConsistent = True # False
    if UseSelfConsistent:
        T_wt = MSMLib.estimate_transition_matrix(C_wt)  # row-normalized estimate
        T_mut = MSMLib.estimate_transition_matrix(C_mut)  # row-normalized estimate
        T_comb = MSMLib.estimate_transition_matrix(C_comb)  # row-normalized estimate
        T_wt, pops_wt = self_consistent_estimation(T_wt)
        T_mut, pops_mut  = self_consistent_estimation(T_mut)
        T_comb, pops_comb  = self_consistent_estimation(T_comb)
    else:  # use MLE
        T_wt, pops_wt = get_populations(C_wt, method='mle', initial_guess=C_wt)
        T_mut, pops_mut = get_populations(C_mut, method='mle', initial_guess=C_mut)
        T_comb, pops_comb = get_populations(C_comb, method='mle', initial_guess=C_comb)

    # Calculate metrics for ach stateso we can choose the next state to sample 
    JSD, var_JSD = [], []
    surprisals, var_surprisals = [], []
    print 'rows with surprisal == 0:', 
    for k in range(nstates):
     
        mean_surprisal = surprisal(C_wt[k,:], C_mut[k,:], normalize=True)
        # VAV debug
        if mean_surprisal == 0.0:
            print k,  #, 'C_wt[k,:], C_mut[k,:]', C_wt[k,:], C_mut[k,:]

        UseBootstrap = False
        if UseBootstrap:
            var_surprisal  = surprisal_var(C_wt[k,:], C_mut[k,:], normalize=True, bootstrap=True, n_bootstraps=100)
        else:
            var_surprisal  = surprisal_var(C_wt[k,:], C_mut[k,:], normalize=True, bootstrap=False, n_bootstraps=100)

        ### surprisal_var(c1, c2)
        surprisals.append( mean_surprisal )
        var_surprisals.append( var_surprisal )

        # compute the JSD, which we use as a measure to compare the two MSMs

        # approximate JSD using pi_si! 8/2014  it's fast but will it work?
        mean_JSD_k = pops_comb[k]*mean_surprisal
        var_JSD_k = (pops_comb[k]**2)*var_surprisal   # added square!!!

        #mean_JSD_k, var_JSD_k = estimate_JSD(C_wt[k,:], C_mut[k,:], pops_wt[k], pops_mut[k], pops_comb[k])
        #JSD_k = get_JSD(C_wt[k,:], C_mut[k,:], pops_wt[k], pops_mut[k], pops_comb[k])
        JSD.append(mean_JSD_k)
        var_JSD.append(var_JSD_k)

    sum_JSD, var_sum_JSD = np.array(JSD).sum() , np.array(var_JSD).sum()
    sum_surprisal, var_sum_surprisal = np.array(surprisals).sum() , np.array(var_surprisals).sum()
    total_JSD.append(sum_JSD)
    total_JSD_var.append(var_sum_JSD)

    #======== output ==========# 
    print '#============================================#'
    print '#state_WT\tstate_MUT\tnsamples\tJSD\tvar(JSD)\tJSD_i...'

    # write jsd data to file
    jsd_fields = [str(j_wt), str(j_mut), str(total_samples), '%e'%sum_JSD, '%e'%var_sum_JSD] + [ '%e'%value  for value in JSD ] + [ '%e'%value  for value in var_JSD ]
    jsd_line = string.joinfields(jsd_fields, '\t')
    fout_jsds.write(jsd_line+'\n')
    print jsd_line

    # write surprisal data to file
    surprisal_fields = [str(j_wt), str(j_mut), str(total_samples), '%e'%sum_surprisal, '%e'%var_sum_surprisal] + [ '%e'%value  for value in surprisals ] + [ '%e'%value  for value in var_surprisals ]
    surprisal_line = string.joinfields(surprisal_fields, '\t')
    fout_surprisals.write(surprisal_line+'\n')
    print 's_i:\n\t', surprisal_line

    # write pops data to file
    pops_fields = [str(j_wt), str(j_mut), str(total_samples)] + [ '%e'%value  for value in pops_comb ] 
    pops_line = string.joinfields(pops_fields, '\t')
    fout_pops.write(pops_line+'\n')
    print 'pops:\n\t', pops_line

    if (adaptive_method == 'surprisal'):
        # Pick the largest surprisal as the next state to sample
        j_wt = np.argmax( surprisals )
        j_mut = j_wt 

    elif (adaptive_method == 'surprisal variance'):
        # Pick the largest surprisal as the next state to sample
        j_wt = np.argmax( var_surprisals )
        j_mut = j_wt

    elif (adaptive_method == 'JSD'):
        # Pick the largest JSD as the next state to sample
        j_wt = np.argmax( JSD )
        j_mut = j_wt
        print 'picked largest JSD =', j_wt

    elif adaptive_method == 'JSD variance':
        # Pick the largest var(JSD) i.e. uncertainty as the next state to sample
        j_wt = np.argmax( var_JSD )
        j_mut = j_wt

    elif adaptive_method == 'random':
        # Pick the next state at random
        j_wt = np.random.randint(nstates)
        j_mut = j_wt

    elif adaptive_method == 'continue traj': 
        # Continue the trajectory from where either wt or mut left off
        j_wt = final_state_wt
        print 'continuing from final wt state', j_wt
        j_mut = final_state_mut
        print 'continuing from final mut state', j_mut
        
    elif adaptive_method == 'pi_si':
        j_wt = np.argmax( pops_comb*surprisals )
        j_mut = j_wt

    elif adaptive_method == 'pi_si variance':
        j_wt = np.argmax( pops_comb*var_surprisals )
        j_mut = j_wt

    elif adaptive_method == 'nina':  # Nina Singhal's method of sensitivity analysis (SA)
        ssa = ssaCalculator(1, C_comb)  # use combined counts for the analysis
        #""" Class for objects to govern sensitivity calculations. """
        #def __init__( self, lagTime, tcMatrix, priorCounts=.5, evalList=[1], nNewSamples=1, timeUnits="ns", recommendationScheme = 'Nina' )
        eigIndex = 0
        resample0 = ssa.resamplingDistribution(0)  # returns [ 0 0 0 .... 0 1 0 ... 0 ] for sampled index
        j_wt = np.argmax(np.array(resample0))
        j_mut = j_wt


# Close open file handles
fout_surprisals.close()
fout_jsds.close()
fout_pops.close()


#if UsePyplot:
if (1):
    # plot the results
    plt.figure()
    time = np.arange(0,total_samples, nsamples)

    plt.subplot(1,2,1)    
    plt.errorbar(time, total_JSD, yerr=np.array(total_JSD_var)**0.5)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('time (steps)')
    plt.ylabel('total_JSD')

    plt.subplot(1,2,2)
    plt.loglog(time, total_JSD_var, 'k-')
    plt.xlabel('time (steps)')
    plt.ylabel('total_JSD_var')

    print 'Made it here!'

    plt.savefig(outname+'.pdf')
    #print 'Wrote', outname+'.pdf', '...'
    #plt.show()

