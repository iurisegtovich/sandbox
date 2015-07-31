#! /usr/bin/env python

import random
import string
import math
import os, sys, copy, pickle

import scipy.sparse
import scipy.linalg
import scipy

from msmbuilder.msm_analysis import *
from msmbuilder.MSMLib import *

from scipy.io import mmread, mmwrite


usage = """Usage: reweightT.py infile.mtx outfile_micro.mtx outfile_macro.mtx microstates.dat sequence 

    sequence - a one-letter amino acid sequence, for the six hydrophobic(H) residues
               in the HPHPHPHPPHPH sequence:  AAAAAA  for example

    Will write files [sequence].dat containing the stability and slowest, next-slowest implied timescales
    """


sys.path.append('../../scripts')
from HelperTools  import *


# Functions

def get_timescales_and_eigs(tprob, num_modes = 10, lag_time = 1., Verbose=False):
    if Verbose: print 'Computing Implied Timescales....'

    if num_modes + 1 > tprob.shape[0]:
        print "cannot get %d eigenmodes from a rank %d matrix"%(num_modes + 1, tprob.shape[0])
        print "Getting as many modes as possible..."
        #logger.warning("cannot get %d eigenmodes from a rank %d matrix", num_modes + 1, tprob.shape[0])
        #logger.warning("Getting as many modes as possible...")
        num_modes = tprob.shape[0] - 1
    
    eigenvalues, eigenvectors = get_eigenvectors(tprob, num_modes + 1, right=False, tol=1e-30)

    # discard the stationary eigenmode
    teigenvalues = np.real(eigenvalues[1:])
    teigenvectors = np.real(eigenvectors[:, 1:])
    
    implied_timescales = - lag_time / np.log(teigenvalues)
    lagtimes = lag_time * np.ones((num_modes)) 

    result =  np.zeros((num_modes,2))
    result[:,0] = lagtimes
    result[:,1] = implied_timescales

    return result, eigenvalues, eigenvectors  


def reweightTransitionMatrix(T, sequence, MJdict, ContactStateIndices, microContactStates, beta, Verbose=False):
    if Verbose: print 'reweighting the transition matrix to reflect the sequence...'
    newT = T.tolil()
    rows, cols = newT.nonzero()
    for k in range(len(rows)):
        if k%10000 == 0:
           if Verbose:  print k, 'of', len(rows), 'nonzero entries'
        i, j = rows[k], cols[k]
        if ContactStateIndices[i] != ContactStateIndices[j]:
            u_i = energy(microContactStates[i], sequence, MJdict)
            u_j = energy(microContactStates[j], sequence, MJdict)
            scalefactor = np.exp(beta*(u_i-u_j))
            newT[i,j] = newT[i,j]*min(1., scalefactor)

    print 'Estimating equil pops from stat mech...'
    equilpops = np.zeros((len(microContactStates)))
    for k in range(len(microContactStates)):
        equilpops[k] = np.exp(-beta*energy(microContactStates[k], sequence, MJdict))
    equilpops = equilpops/equilpops.sum()

    if Verbose:
        for k in range(len(microContactStates)):
            print k, equilpops[k]    

    return estimate_transition_matrix(newT).tolil(), equilpops





################
# Main program

if len(sys.argv) < 6:
    print usage
    sys.exit(1)

microTFn            = sys.argv[1]
reweighted_microTFn = sys.argv[2]
reweighted_macroTFn = sys.argv[3]
microstatesFn       = sys.argv[4]
sequence            = sys.argv[5]
outfile = sequence + '.dat'

# Read in transition matrix
if os.path.exists(microTFn):
    T = mmread(microTFn)
    print T
else:
    print "Can't find file:", microTFn, '...exiting'
    sys.exit(1)

# Read in microstate info
microstates = []         # list of vecs
microContactStates = [] # list of contact states for each 
microNumContacts  = []   # list of number of contacts

m =  MicrostateInfo(microstatesFn)
print 'native index:', m.NativeContactStateIndex


# Build sequence information
#if os.path.exists(reweighted_microTFn):
#        newT = mmread(reweighted_microTFn)
#        print newT
#else:
if (1):
        MJdict = build_MJ_dict()
        residues = "Cys Met Phe Ile Leu Val Trp Tyr Ala Gly Thr Ser Asn Gln Asp Glu His Arg Lys Pro".split()
        resIndices = [0,2,4,6,9,11]
        seqdict = {} # {resindex: 'Cys'} e.g.
        initial_sequence = [olc2tlc(sequence[i]) for i in range(len(sequence))]
        for i in resIndices:
            seqdict[i] = initial_sequence[resIndices.index(i)]
        print 'seqdict = ', seqdict
        beta = 1.0  # because MJ units are in kT
        newT, new_equilpops_from_statmech = reweightTransitionMatrix(T, seqdict, MJdict, m.ContactStateIndices, m.microContactStates, beta)
        mmwrite(reweighted_microTFn, newT, comment=sequence)

CalculateOriginalTimescales = True #False

if CalculateOriginalTimescales:
    print 'Original T:'
    result, evals, evecs = get_timescales_and_eigs(T, num_modes = 10)
    print evals, evecs
    print result

    # Test to see if the equilpops match our stat mech reasoning 
    equilpops = evecs[:,0]/evecs[:,0].sum()
    epsilon = 0.0  # in units kT
    equilpops_from_theory = np.exp(-epsilon*np.array(m.microNumContacts))
    equilpops_from_theory = equilpops_from_theory/equilpops_from_theory.sum()
    print 'equilpops', equilpops
    print 'equilpops_from_theory', equilpops_from_theory

    # calculate the stability of the native state
    INative = np.array(m.ContactStateIndices)==m.NativeContactStateIndex
    stability = np.real( evecs[INative,0]/evecs[:,0].sum() )
    # print stability and timescale info 
    print '#stability\tslowest\tnext-slowest'
    print '%8.3f\t%16.8f\t%16.8f'%(stability, result[0,1], result[1,1])


print 'reweighted T:'
result, evals, evecs = get_timescales_and_eigs(newT, num_modes = 10)
print evals, evecs
print result

if (1):
    # Test to see if the equilpops match our stat mech reasoning 
    equilpops = evecs[:,0]/evecs[:,0].sum()
    equilpops_from_theory = new_equilpops_from_statmech
    print 'equilpops', equilpops
    print equilpops[11805:11820]
    print 'equilpops_from_theory', equilpops_from_theory
    from matplotlib import pyplot as plt
    plt.figure()
    plt.subplot(2,2,1)
    plt.plot(equilpops, equilpops_from_theory,'.')
    plt.title('microstate populations')
    plt.xlabel('eqpops from MSM')
    plt.ylabel('eqpops from Boltzmann')
    plt.xscale('log')
    plt.yscale('log')


# calculate the stability of the native state
INative = np.array(m.ContactStateIndices)==m.NativeContactStateIndex
stability = np.real( evecs[INative,0]/evecs[:,0].sum() )

# print stability and timescale info 
header = '#sequence\tstability\tslowest\tnext-slowest'
line = '%s\t%8.3f\t%16.8f\t%16.8f'%(sequence, stability, result[0,1], result[1,1])
print header+'\n'+line

# write stability and timescale info to output file
fout = open(outfile,'w')
fout.write(header+'\n')
fout.write(line+'\n')
fout.close()


############################################################
### Build a macrostate version of the transition matrix  ###

# Calculate the timescale and stability of the macrostate transition matrix

NumContactStates = len(m.uniqueContactStates)
C_macro  = scipy.sparse.lil_matrix((int(NumContactStates),int(NumContactStates)))

newT_micro = newT.tolil()
print 'newT_micro.rows', newT_micro.rows
print 'newT_micro.data', newT_micro.data

for i in range(newT.shape[0]):
    for j in range(len(newT_micro.rows[i])):
        k = newT_micro.rows[i][j]
        C_macro[m.ContactStateIndices[i], m.ContactStateIndices[k]] += equilpops_from_theory[i]*newT_micro.data[i][j]

T_macro = estimate_transition_matrix(C_macro) # this gives the naive answer (not MLE)
print 'T_macro', T_macro

# write new macro to file
mmwrite(reweighted_macroTFn, T_macro, comment=sequence)


print 'reweighted T_macro:'
result, evals, evecs = get_timescales_and_eigs(T_macro, num_modes = 10)
print evals, evecs
print result

# calculate the stability of the native state from the MACROSTATE 
INative = np.array(range(int(NumContactStates)))==m.NativeContactStateIndex
stability = np.real( evecs[INative,0]/evecs[:,0].sum() )
# print stability and timescale info 
print '#stability\tslowest\tnext-slowest'
print '%8.3f\t%16.8f\t%16.8f'%(stability, result[0,1], result[1,1])

if (1):
    # Test to see if the equilpops match our stat mech reasoning 
    macro_equilpops = evecs[:,0]/evecs[:,0].sum()

    macro_equilpops_from_theory = np.zeros(int(NumContactStates))
    for k in range(len(m.ContactStateIndices)):
        macro_equilpops_from_theory[m.ContactStateIndices[k]] += equilpops_from_theory[k] 

    print 'macro_equilpops_from_theory.sum()', macro_equilpops_from_theory.sum()

    print 'macro_equilpops', macro_equilpops
    print 'macro_equilpops_from_theory', macro_equilpops_from_theory

    I1 = np.argsort(-macro_equilpops)
    #I2 = np.argsort(-macro_equilpops_from_theory)
    print '# state\tmacro_equilpops\tmacro_equilpops_from_theory'
    for k in range(10):
        print I1[k], macro_equilpops[I1][k], macro_equilpops_from_theory[I1][k]
     
    plt.subplot(2,2,2,)
    plt.plot(macro_equilpops, macro_equilpops_from_theory,'.')
    plt.title('macrostate populations')
    plt.xlabel('eqpops from MSM')
    plt.ylabel('eqpops from Boltzmann')
    plt.xscale('log')
    plt.yscale('log')
    plt.show()




