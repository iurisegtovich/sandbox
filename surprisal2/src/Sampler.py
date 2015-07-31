# New Sampler object Aug 2014
# GOAL is to have the same object to sample sparse and dnse matrices! 

import os, sys, glob, string, copy, pickle
import math, random
import pickle


import scipy.sparse
import scipy.linalg
import scipy
from scipy.io import mmread, mmwrite
import numpy as np

try:
    from msmbuilder.MSMLib import ergodic_trim
except:
    from msmbuilder.MSMLib import ErgodicTrim as ergodic_trim

""" FIX THIS -- This shouls be a method or application of the sampler object --VAV 8.2014
def ensure_ergodicity(sampler_obj, count_mtx, num_samples=25):
    trimmed_counts, Mapping = ergodic_trim(scipy.sparse.csr_matrix(count_mtx)) 
    ergodic_mtx = False
    i = 0        
    while not ergodic_mtx:
        i += 1
        if trimmed_counts.shape == count_mtx.shape:
            ergodic_mtx = True
            continue
        if i % 1000 == 0:
            print "Attempt %d for Ergodicity"%(i)
        nonergodic_index = (Mapping == -1)
        sample_dist = np.where(nonergodic_index, num_samples, 0) 
        count_mtx = sampler_obj.sample_from_distribution(sample_dist, count_mtx) 
        trimmed_counts, Mapping = ergodic_trim(scipy.sparse.csr_matrix(count_mtx)) 
    return count_mtx

"""

def draw_index(probs, n_picks=1, UseFastMethod=True):
    """Draw a number (or many numbers, controlled by n_picks), weighted by the probabilities probs."""

    if UseFastMethod:

        t = np.cumsum(probs)
        s = sum(probs)
        return np.searchsorted(t,np.random.rand(n_picks)*s)

    else:
      try:
        nprobs = len(probs)
      except:
        nprobs = probs.shape[0]
      r = np.random.random()
      i = 0
      pthresh = 0.0
      while i < nprobs:
        pthresh += probs[i]
        if r < pthresh:
           return i
        i += 1
      return min(nprobs-1, i)




class Sampler:

    def __init__(self, tProb, jump_Fn):

        # make sparse matrices LIL
        if scipy.sparse.issparse(tProb):
            self.tProb = tProb.tolil()
            self.issparse = True
            self.nstates = tProb.shape[0]
        # otherwise store a dense matrix
        else:
            self.tProb = tProb
            self.issparse = False
       
        # Build or read jump probabilities      
	if os.path.exists(jump_Fn):
            print 'Reading jump probabilities from: ', jump_Fn, '...'
            fin = open(jump_Fn, 'r')
            self.jumps = pickle.load(fin)
	    fin.close()
            print '... Done!'
    	else:
            self.jumps = self.find_jumps(self.tProb, jump_Fn)

    def find_jumps(self, tProb, jump_Fn):
        """Creates a dictionary of jump probabilities. Every elements of
           jumps is the (possible transition, probability of transition)"""

        nstates = tProb.shape[0]
        print 'Making a dictionary of jump probabilities...'
        jumps = {}

        if self.issparse: 
            for i in range(nstates):
                possible_j = np.array(tProb[i, :].nonzero()[1], dtype=np.int)
                j_probs = np.array(tProb[i, possible_j].todense())
                j_probs = j_probs.reshape(j_probs.shape[1],)
                jumps[i] = (possible_j, j_probs)
                if i%100 == 0:
                    print i, jumps[i]
        else:
            for i in range(nstates):
                possible_j = np.array(self.tProb[i, :].nonzero()[0], dtype=np.int)
                j_probs = self.tProb[i, possible_j]
                jumps[i] = (possible_j, j_probs)
                if i%100 == 0:
                    print i, jumps[i]

        print 'Pickling jumps to ', jump_Fn, '...'
        fout = open(jump_Fn, 'w')
        pickle.dump(jumps, fout)
        fout.close()
        return jumps



        print 'Pickling jumps to ', jump_Fn, '...'
        fout = open(jump_Fn, 'w')
        pickle.dump(jumps, fout)
        fout.close()
        return jumps


    def sample(self, state, nsamples, previous_counts=None):
        """Samples transition counts from a particular state

        ARGUMENTS
          state - the index of the state to sample from
          nsamples  - the number of transition countis to be samples from that state

        OPTIONAL ARGUMENTS
          previous_counts - if supplied, counts will be added to previously-supplied count matrix 

        RETURNS
          counts - a dense (!) 2D array of counts 

        """

        if previous_counts == None:
            counts = np.zeros( self.tProb.shape , dtype=np.float64)
        else:
            # check to see if previous_counts is the right shape
            np.testing.assert_array_equal(self.tProb.shape, previous_counts.shape, err_msg='The previously supplied counts are not the right shape', verbose=True)
            counts = np.copy(previous_counts)

        # draw counts
        counts[state,:] += np.random.multinomial(nsamples, self.tProb[state,:])

        return counts



    def sample_uniform(self, num_samples, new_count_mtx=None):
        """Does uniform sampling by calling the sample_from_distribution
           method"""
        num_states = self.tProb.shape[0]
        sample_dist = [num_samples for i in range(num_states)]

        if new_count_mtx == None:
            new_count_mtx = self.sample_from_distribution(sample_dist)
        else:
            new_count_mtx = self.sample_from_distribution(sample_dist, new_count_mtx)
        return new_count_mtx



    def sample_traj(self, state, nsamples, previous_counts=None):
        """Samples counts from a Markov chain trajectory starting from a particular state

        ARGUMENTS
          state - the index of the state to start sampling from
          nsamples  - the number of transition countis to be samples from that state

        OPTIONAL ARGUMENTS
          previous_counts - if supplied, counts will be added to previously-supplied count matrix 

        RETURNS
          counts - a new array of counts 
          final_state - the last state in the trajectory to be sampled (so we can choose to continue from that point)

        """
        
        # set up sampling for sparse arrays 
        if self.issparse:
            num_states = self.tProb.shape[0]
            if previous_counts == None:
                counts = scipy.sparse.lil_matrix((int(num_states), int(num_states)))
            else:
                counts = copy(previous_counts)
   
        # set up sampling for dense arrays    
        else:
            counts = np.zeros( self.tProb.shape , dtype=np.float64)
            if not (previous_counts == None):
                # check to see if previous_counts is the right shape
                np.testing.assert_array_equal(self.tProb.shape, previous_counts.shape, err_msg='The previously supplied counts are not the right shape', verbose=True)
                counts += previous_counts

        # do the sampling
        for i in range(nsamples):
            newstate = self.weighted_pick(self.tProb[state,:], 1)[0]
            counts[state,newstate] += 1.0
            state = newstate

        return counts, newstate


    def weighted_pick(self, weights, n_picks):
        """
        Weighted random selection returns n_picks random indexes.  the chance to pick the index i is give by the weight weights[i].
        """

        t = np.cumsum(weights)
        s = sum(weights)
        return np.searchsorted(t,np.random.rand(n_picks)*s)


    def sample_from_distribution(self, sample_dist, new_count_mtx=None):
        """Returns  a count matrix given some sampling distribution based on the 
           prevously found jump probabilities. If the third arguement is given
           it will add to a previous count matrix rather than create a new one."""

        num_states = self.tProb.shape[0]

        if (len(sample_dist) != num_states):
            print 'Distribution does not correspond to probabilities. Try again.'
            sys.exit(1)

        if new_count_mtx == None:
            new_count_mtx = scipy.sparse.lil_matrix((int(num_states), int(num_states)))
        else:
            new_count_mtx = new_count_mtx.tolil()

        total_counts = 0

        for i in np.array(sample_dist).nonzero()[0]:
            #if i%1000 == 0:
                #print 'Sampling ', sample_dist[i], 'transitions from state ', i, ' of ', num_states, '...'
            possible_j = self.jumps[i][0]
            j_probs = self.jumps[i][1]
            picks = draw_index(j_probs, n_picks=sample_dist[i])
            picked_counts = np.bincount(picks, minlength=len(possible_j))
            for k in range(picked_counts.shape[0]):
                new_count_mtx[i, possible_j[k]] += picked_counts[k]

            total_counts += sum(sample_dist)
            #print 'Sampled a total of %d counts from %d states'%(total_counts, num_states)

        return new_count_mtx


