#!/usr/bin/env python
"""Last updated: 01/03/2012. Sampler class for the creation of new 
   count matrices from a prior transition probability matrix. """
import random
import string
import math
import os
import sys
import copy 
import pickle
import scipy.sparse
import scipy.linalg
import scipy
from scipy.io import mmread, mmwrite
import numpy as np
from HelperTools import draw_index
try:
    from msmbuilder.MSMLib import ergodic_trim
except:
    from msmbuilder.MSMLib import ErgodicTrim as ergodic_trim

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

class Sampler:
    def __init__(self, tProb, jump_Fn):
        self.tProb = tProb.tolil()

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
        num_states = tProb.shape[0]
        print 'Making a dictionary of jump probabilities...'
        jumps = {}

        for i in range(num_states):
            possible_j = np.array(tProb[i, :].nonzero()[1], dtype=np.int)
            j_probs = np.array(tProb[i, possible_j].todense())
            j_probs = j_probs.reshape(j_probs.shape[1],)
            jumps[i] = (possible_j, j_probs)
            if i%100 == 0:
                print i, jumps[i]

        print 'Pickling jumps to ', jump_Fn, '...'
        fout = open(jump_Fn, 'w')
        pickle.dump(jumps, fout)
        fout.close()
        return jumps

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

    def sample_traj(self, num_samples, initial_state=0, count_mtx=None):
        state = initial_state
        num_states = self.tProb.shape[0]

        if count_mtx == None:
            count_mtx = scipy.sparse.lil_matrix((int(num_states), int(num_states)))
        else:
            count_mtx = count_mtx.tolil()
    
        for i in range(num_samples):
            possible_j = self.jumps[state][0]
            jprobs = self.jumps[state][1]
            new_state = possible_j[draw_index(jprobs, n_picks=1)[0]]
            count_mtx[state, new_state] += 1
            state = new_state
        return count_mtx, state
 




if __name__ == '__main__':
    usage = """Sampler.py transition_prob_mtx jump_File_name num_sample. This will do uniform sampling based on transition_prob_mtx for num_samples for each state. Jumps are read or written to the file named jump_File_name"""
    tProb = mmread(sys.argv[1])
    jumpFn = sys.argv[2]

    if len(sys.argv) != 4:
        print usage

    model = Sampler(tProb, jumpFn)
    new_count_mtx = model.sample_uniform(int(sys.argv[3]))
    mmwrite("tCounts.%d.mtx"%(num_samples), new_count_mtx) 
