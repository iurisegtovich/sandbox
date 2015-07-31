#/usr/bin/env python

"""
SurprisalCalculator objects are used to calculate the surprisal
values for two sparse transition count matrices. The variance of 
the surprisal values can also be calculated.
"""

import sys
import os 
import numpy as np
from scipy.io import mmread, mmwrite
from scipy.sparse import isspmatrix_coo	
from sklearn.preprocessing import normalize
from msmbuilder.MSMLib import estimate_transition_matrix, mle_reversible_count_matrix  #, __mle_reversible_count_matrix_lutz__ #For population calculation
from msmbuilder import msm_analysis
from Reversible_MLE_with_Guess import *

def get_populations(count_mtx,  method='mle', initial_guess=None, Verbose=False):
    """ Get the populations (first eigenvector of the transition matrix). Necessary
        for JSD variance calculations. There are two possible symmetrization methods 
        to estimate transition matrices from count matrices, namely MLE and a naive
        count-based symmetrization. """

    # Uses naive symmetrization method to determine the populations. 
    # Do not use exclusively if the count matrix does not already reflect
    # the equilibrium populations. 
    if method == 'counts':
        counts = count_mtx.todense()
        counts_sym = 0.5*(counts.transpose()+counts)
        return np.array(counts_sym.sum(axis=1)/counts_sym.sum()).flatten()

    # Uses Maximum Likelihood Estimator from MSMBuilder2 to find the symmetrized
    # counts matrix. Initial guess should be the symmetrized matrix from the 
    # previous step. 
    if method == 'mle':
        if (initial_guess == None): 
            counts_sym = mle_reversible_count_matrix_with_guess(count_mtx)
        else:
            counts_sym = mle_reversible_count_matrix_with_guess(count_mtx, 
                                                                initial_guess)
        T = estimate_transition_matrix(counts_sym)
        evals, evecs = msm_analysis.get_eigenvectors(T, 6, 0.001, 50)
        return np.real(evecs[:,0]/evecs[:,0].sum()), counts_sym

def get_harmonic_mean_population(count_mtx1, count_mtx2):
    try:
        populations1, self.new_init_guess1 = get_populations(count_mtx1, 'mle', count_mtx1)
        populations2, self.new_init_guess2 = get_populations(count_mtx2, 'mle',  count_mtx2)
    except:  
        populations1 = get_populations(count_mtx1, 'counts')
        populations2 = get_populations(count_mtx2, 'counts')
    return (2*populations1*populations2)/(populations1+populations2) 

class SurprisalCalculator:
    """ SurprisalCalculator objects store information concerning the current 
        populations and surprisal values. Requires two counts matrices as input
        along with either previous symmetrized count matrices (if pop_method is
        'mle') or previous populations (if pop_method is 'retain')  """ 
    def __init__(self, sparse1, sparse2, pop_method=None, init_guess1=None, init_guess2=None, init_guess_comb=None):

        self.sparse1 = sparse1.tocoo()
        self.sparse2 = sparse2.tocoo()
        self.comb_matrix = (self.sparse1 + self.sparse2).tocoo()

        if pop_method == None:
            pop_method = 'counts'   # Set the default population method
   
        if(pop_method == 'retain_first_mle'):
            if(init_guess2 == None):
                self.populations2, self.new_init_guess2 = get_populations(self.sparse2, 
                                                                          'mle', self.sparse2)
            else:
                self.populations2, self.new_init_guess2 = get_populations(self.sparse2, 
                                                                          'mle', init_guess2)
            self.populations1 = init_guess1

        if (pop_method == 'retain'):
            self.populations1 = init_guess1 #"init_guesses" should be previous 
                                            #populations for retain!
            self.populations2 = init_guess2
            self.comb_pops = init_guess_comb

        if (pop_method == 'counts'):
            self.populations1 = get_populations(self.sparse1, 'counts')
            self.populations2 = get_populations(self.sparse2, 'counts')
            self.comb_pops = get_populations(self.sparse1+self.sparse2, 'counts') 

        if (pop_method == 'mle'): 
            if(init_guess_comb == None):
                self.comb_pops, self.new_init_guess_comb = get_populations(self.comb_matrix, 
                                                                           'mle', self.comb_matrix)
            else:
                self.comb_pops, self.new_init_guess_comb = get_populations(self.comb_matrix, 
                                                                           'mle', init_guess_comb)
            if(init_guess1 == None):
                self.populations1, self.new_init_guess1 = get_populations(self.sparse1, 
                                                                          'mle', self.sparse1)
            else:
                self.populations1, self.new_init_guess1 = get_populations(self.sparse1, 
                                                                          'mle', init_guess1)
            if(init_guess2 == None):
                self.populations2, self.new_init_guess2 = get_populations(self.sparse2, 
                                                                          'mle', self.sparse2)
            else:
                self.populations2, self.new_init_guess2 = get_populations(self.sparse2, 'mle', init_guess2)


        self.exceptions = 0

        if (self.sparse1.shape != self.sparse2.shape):
            print "ERROR: Systems much have the same discretized states." 
            print "       The matrices must be of the same shape."
            sys.exit(1)
 
    def prepare_count_arrays(self, i):

        """
        Prepares count arrays from sparse matrices within a  
        SurprisalCalculator object for surprisal calculations
        for every state i. 

        Returns two arrays of counts with 0 inserted for states
        which aren't shared between the models. 
        """
        sparse1 = self.sparse1
        sparse2 = self.sparse2

        j1 = sparse1.col[sparse1.row == i] 
        j2 = sparse2.col[sparse2.row == i] 
        j_states = np.unique(np.append(j1, j2))
        j_states.sort()

        counts1 = sparse1.data[sparse1.row == i]
        counts2 = sparse2.data[sparse2.row == i] 

        for index in range(len(j_states)):
            if not j1.__contains__(j_states[index]):  
                counts1 = np.insert(counts1, index, 0)     
            if not j2.__contains__(j_states[index]):  
                counts2 = np.insert(counts2, index, 0)

        return counts1.astype(np.float64), counts2.astype(np.float64)

    def calculate_entropy(self, counts):
        """Calculates the entropy of a given set of transition counts, which is
           S_i = sum_j n_{ij}/N_i*ln(n_{ij}/N_i)"""
        entropy = 0.0
        total_counts = np.sum(counts).astype(np.float64)
        for index in range(len(counts)):
            if counts[index] != 0:
                entropy -= (counts[index]/total_counts)*np.log(counts[index]/total_counts)
        if total_counts == 0:
            entropy = 0.0
        return entropy

    def calculate_surprisal(self, counts1, counts2, normalize=False):
        """
        Calculates surprisal for two sets of counts. If normalized, the 
        surprisal value is divided by the total number of counts from that state. 
        """
  
        combined_counts = counts1 + counts2
        total_combined_counts = np.sum(combined_counts)
        total_counts1 = np.sum(counts1)
        total_counts2 = np.sum(counts2)

        combined_entropy = self.calculate_entropy(combined_counts)
        entropy1 = self.calculate_entropy(counts1)
        entropy2 = self.calculate_entropy(counts2)

        surprisal = ((total_combined_counts * combined_entropy) -
                        (total_counts2 * entropy2) - (total_counts1 * entropy1))

        if surprisal < 0:
            print "Warning! Negative Surprisal value! Defaulting to 0!"
            surprisal = 0.0

        if (normalize and total_combined_counts != 0):
            surprisal /= float(total_combined_counts)

        return surprisal
 
    def calculate_all_surprisal(self, verbose=False, normalized=False, weighted=False):
        """
        Calculates surprisal for all states of the sparse matrices
        in the SurprisalCalculator object. Returns an array where
        surprisals[i] is the surprisal of state i. 
    
        The surprisal values can be both normalized by the number 
        of counts and weighted by the populations.
        """
      
        surprisals = []
        for i in range(self.sparse1.shape[0]):
            if (verbose):
                print("Working on state %d"%i) 

            counts1, counts2 = self.prepare_count_arrays(i) 

            if(not normalized):
                surprisals.append(self.calculate_surprisal(counts1, counts2))
            if(normalized):
                surprisals.append(self.calculate_surprisal(counts1, counts2, True))

        if (weighted):
            try:
                populations = get_harmonic_mean_population(self.sparse1, self.sparse2)
                surprisals *= populations 
            except:
                self.exceptions += 1
                mmwrite("Failed_WT%d.mtx"%(self.exceptions), self.sparse1)
                mmwrite("Failed_MUT%d.mtx"%(self.exceptions), self.sparse2)
                if (self.exceptions == 5):
                    print "Five exceptions. Closing program"
                populations = get_populations(self.sparse2, 'counts')
                
        return surprisals

    def get_covariance_matrix(self, sparse):
        """
            Returns the covariance matrix of a multinomial distribution for 
            p_i ~ counts
            """
        covariance_matrix = np.zeros((sparse.shape[0], sparse.shape[1]))
        
        for i in range(sparse.shape[0]):
            j_states = sparse.col[sparse.row == i]
            counts = sparse.data[sparse.row == i].astype(np.float64)
            total_counts = np.sum(counts).astype(np.float64) 
            self_transition_counts = sparse.data[(sparse.row == i)*(sparse.col == i)]
            
            if not self_transition_counts: #breaks this iteration of the loop
                continue                   #if there are 0 self transition counts
            
            else:
                self_transition_counts = self_transition_counts[0] 
            
            for index in range(len(j_states)):
                j_state = j_states[index]
                covariance_matrix[i][j_state] = -self_transition_counts*counts[index]/total_counts
                if i == j_state:
                    covariance_matrix[i][i] += self_transition_counts
        return covariance_matrix

    def get_JSD_impacts(self, sparse1, sparse2, nsamples1, nsamples2):
        """Gives predicted changes in JSD_i if new counts nsamples1, nsamples2 are to be sampled from each state i ."""

        if not isspmatrix_coo(sparse1):
            sparse1 = sparse1.tocoo()
        if not isspmatrix_coo(sparse2):
            sparse2 = sparse2.tocoo()

        dHcomb_dn1, dH1_dn1, dHcomb_dn2, dH2_dn2 = self.get_entropy_derivatives(sparse1, sparse2)

        f1 = np.zeros( (sparse1.shape[0], sparse1.shape[1]) ) # a matrix of normalized rows
        for i in range(sparse1.shape[0]):
            # get the normalized row i as from the COO matrix (yeesh!)
            j_states1 = sparse1.col[sparse1.row == i]
            counts1_ij = sparse1.data[sparse1.row == i].astype(np.float64)
            total_counts1 = np.sum(counts1_ij).astype(np.float64)
            for index in range(len(j_states1)):
                j = j_states1[index]
                f1[i][j] =  counts1_ij[index]/total_counts1

        f2 = np.zeros( (sparse2.shape[0], sparse2.shape[1]) ) # a matrix of normalized rows
        for i in range(sparse2.shape[0]):
            # get the normalized row i as from the COO matrix (yeesh!)
            j_states2 = sparse2.col[sparse2.row == i]
            counts2_ij = sparse2.data[sparse2.row == i].astype(np.float64)
            total_counts2 = np.sum(counts2_ij).astype(np.float64)
            for index in range(len(j_states2)):
                j = j_states2[index]
                f2[i][j] =  counts2_ij[index]/total_counts2

        print 'f1', f1
        print 'f2', f2 

        # The expected value of the change in JSD_i for new counts M1 and M2  is
        #
        #      M1 < d JSD_i / dn1_ij > + M2 < d JSD_i / dn2_ij >
        #
        # =    M1 * (pi_comb_i <dHcomb_dn1>_i - pi1_i * (M1/(M1+M2)) <dH1_dn1>_i )
        #    + M2 * (pi_comb_i <dHcomb_dn2>_i - pi2_i * (M2/(M1+M2)) <dH2_dn2>_i ) 


        if self.combined_populations == None:
            self.combined_populations = (self.populations1 + self.populations2)/2.0
        print 'self.populations1', self.populations1
        print ' self.populations2',  self.populations2
        print 'self.combined_populations', self.combined_populations
        M1, M2 = nsamples1, nsamples2

        print 'np.dot(f1[i,:],dHcomb_dn1[i,:])', np.dot(f1[i,:],dHcomb_dn1[i,:])
        print 'np.dot(f1[i,:],dH1_dn1[i,:])', np.dot(f1[i,:],dH1_dn1[i,:])
        print 'np.dot(f2[i,:],dHcomb_dn2[i,:])', np.dot(f2[i,:],dHcomb_dn2[i,:])
        print ' np.dot(f2[i,:],dH2_dn2[i,:])', np.dot(f2[i,:],dH2_dn2[i,:])

        results = []
        for i in range(sparse1.shape[0]): 
            result =   M1 * (self.combined_populations[i] * np.dot(f1[i,:],dHcomb_dn1[i,:]) - self.populations1[i]*M1/(M1+M2)* np.dot(f1[i,:],dH1_dn1[i,:])) \
                    +  M2 * (self.combined_populations[i] * np.dot(f2[i,:],dHcomb_dn2[i,:]) - self.populations2[i]*M2/(M1+M2)* np.dot(f2[i,:],dH2_dn2[i,:]))
            results.append(result)

        return results

        
    
    def get_entropy_derivatives(self, sparse1, sparse2):
        """Returns the derivatives dH_i/dn_ij for each state i
 
        INPUTS
        sparse1		-- sparse matrix of counts for wt.
        sparse2         -- sparse matrix of counts for mut.

        RETURNS
        d(Hcomb_i/dn_ij1), d(H1_i/dn_ij1), d(Hcomb_2/dn_ij2), d(H2_i/dn_ij2)
        """

        if not isspmatrix_coo(sparse1):
            sparse1 = sparse1.tocoo()
        if not isspmatrix_coo(sparse2):  
            sparse2 = sparse2.tocoo()

        combined = sparse1 + sparse2
        dHcomb_dn1 = np.zeros((sparse1.shape[0], sparse1.shape[1]))
        dH1_dn1    = np.zeros((sparse1.shape[0], sparse1.shape[1]))
        dHcomb_dn2 = np.zeros((sparse2.shape[0], sparse2.shape[1]))
        dH2_dn2    = np.zeros((sparse2.shape[0], sparse2.shape[1]))

        for i in range(sparse1.shape[0]):
            counts1 = sparse1.data[sparse1.row == i].astype(np.float64) 
            counts2 = sparse2.data[sparse2.row == i].astype(np.float64)
            j_states1 = sparse1.col[sparse1.row == i]
            j_states2 = sparse2.col[sparse2.row == i]
            print 'counts1', counts1
            print 'counts2', counts2
            total_counts1 = np.sum(counts1).astype(np.float64) 
            total_counts2 = np.sum(counts2).astype(np.float64)
            total_combined_counts = total_counts1 + total_counts2
            print 'total_counts1', total_counts1, 'total_counts2', total_counts2, 'total_combined_counts', total_combined_counts
            
            for index in range(len(j_states1)):
                j = j_states1[index]
                dH1_dn1[i][j] = ((total_counts1 - counts1[index])/(total_counts1**2))*np.log((total_counts1 - counts1[index])/counts1[index])
                dHcomb_dn1[i][j] = ((total_combined_counts - counts1[index])/total_combined_counts**2)*np.log((total_combined_counts- counts1[index])/counts1[index])

            for index in range(len(j_states2)):
                j = j_states2[index]
                dH2_dn2[i][j] = ((total_counts2 - counts2[index])/total_counts2**2)*np.log((total_counts2 - counts2[index])/counts2[index])
                dHcomb_dn2[i][j] = ((total_combined_counts - counts2[index])/total_combined_counts**2)*np.log((total_combined_counts- counts2[index])/counts2[index])

        print 'dHcomb_dn1', dHcomb_dn1
        print 'dH1_dn1', dH1_dn1
        print 'dHcomb_dn2', dHcomb_dn2
        print  'dH2_dn2',  dH2_dn2

        return dHcomb_dn1, dH1_dn1, dHcomb_dn2, dH2_dn2
    
    def estimate_surprisal_variance_analytical(self, istate, normalized=True):
        """
        Analytical estimation of variance. Still not quite working, use
        bootstrap method instead. 
        """
        covariance1 = self.get_covariance_matrix(self.sparse1)
        covariance2 = self.get_covariance_matrix(self.sparse2)
        size = self.sparse1.shape[1] #Corresponds to number of j-states
                                     
        diagonal_covariances = np.zeros((2*size, 2*size), np.float)
        diagonal_covariances[0:size, 0:size] = covariance1
        diagonal_covariances[size:2*size, size:2*size] = covariance2
        sensitivities1 = np.zeros(size)
        sensitivities2 = np.zeros(size)

        j_states1 = self.sparse1.col[self.sparse1.row == istate]
        j_states2 = self.sparse2.col[self.sparse2.row == istate]
        j_states = np.unique(np.append(j_states1, j_states2))
        counts1, counts2 = self.prepare_count_arrays(istate)
        total_counts1 = np.sum(counts1).astype(np.float64)
        total_counts2 = np.sum(counts2).astype(np.float64)
        total_counts = total_counts1 + total_counts2

        if total_counts == 0:
            variance = 0.0

        for i in range(len(j_states)):
            j_state = j_states[i]

            if counts1[i] != 0:
                sensitivities1[j_state] += (np.log(total_counts1 / total_counts) - 
                                            np.log(counts1[i] / (counts1[i] + counts2[i])))

            if counts2[i] != 0:
                sensitivities2[j_state] += (np.log(total_counts2 / total_counts) -
                                            np.log(counts2[i] / (counts1[i] + counts2[i])))
                
            if normalized and total_counts != 0:
                sensitivities1[j_state] /= total_counts
                sensitivities2[j_state] /= total_counts

            sensitivities = np.append(sensitivities1, sensitivities2)

            variance = (sensitivities.T).dot(diagonal_covariances).dot(sensitivities)

        return variance

    def estimate_surprisal_variance_bootstrap(self, counts1, counts2, n_bootstraps=10, normalized=False):
        """
        Calculates the variance of two sets of transition counts by 
        creating pseudo-count matrices from transition probabilities
        and then estimating the variance of the resulting array.    
        """
 
        total_counts1 = np.sum(counts1)
        if total_counts1 == 0:      #If total counts is zero, 
            prob1 = counts1         #all probabilities are zero, which will be
        else:                       #equivalent to the count matrix.
            prob1 = np.divide(counts1, total_counts1)

        total_counts2 = np.sum(counts2)
        if total_counts2 == 0:              
            prob2 = counts2   
        else:
            prob2 = np.divide(counts2, total_counts2)

        surprisals = []

        #Draw counts from a multinomial distribution based on the
        #transition probabilities 
        sampled_counts1 = np.random.multinomial(total_counts1, prob1,  
                                                size = n_bootstraps)
        sampled_counts2 = np.random.multinomial(total_counts2, prob2,
                                                size = n_bootstraps) 

        for trial in range(n_bootstraps):
            surprisal = self.calculate_surprisal(sampled_counts1[trial,:],
                                                 sampled_counts2[trial,:], 
                                                 normalized)
            surprisals.append(surprisal)

        return np.array(surprisals).var()

    def calculate_all_variance(self, weighted=True, n_bootstraps=50):
        """ Calculates the variance of the surprisal for every state i using the
           bootstrap method""" 
        variances = []
        populations = (1.0/2.0)*(self.populations1+self.populations2)
        for istate in range(self.sparse2.shape[0]):
            counts1, counts2 = self.prepare_count_arrays(istate)
            if weighted == False:
                variance = self.estimate_surprisal_variance_bootstrap(counts1, counts2, n_bootstraps, normalized=True)
            else: 
                variance = populations[istate]*populations[istate]*self.estimate_surprisal_variance_bootstrap(counts1, counts2, n_bootstraps, normalized=True) 
            variances.append(variance)
        return variances

    def calculate_jsd(self, state, counts1, counts2):
        """Calculates Jennsen Shannon Divergence for state with corresponding 
           transition counts from two models. """
        combined_counts = counts1+counts2 
        total_combined_counts = np.sum(combined_counts)
        total_counts1 = np.sum(counts1)
        total_counts2 = np.sum(counts2)

        combined_entropy = self.calculate_entropy(combined_counts)
        entropy1 = self.calculate_entropy(counts1)
        entropy2 = self.calculate_entropy(counts2)
        return (self.comb_pops[state]*combined_entropy)-(1.0/2.0)*(self.populations1[state]*entropy1+self.populations2[state]*entropy2)
         
        if self.combined_populations == None:
            self.combined_populations = (self.populations1 + self.populations2)/2.0
        combined_entropy = self.calculate_entropy(combined_counts)
        entropy1 = self.calculate_entropy(counts1)
        entropy2 = self.calculate_entropy(counts2)

        p1 = total_counts1/(total_counts1+total_counts2)

        # VAV: fixed the definition to be more precise
        return self.combined_populations[state]*combined_entropy - p1*self.populations1[state]*entropy1 - (1.-p1)*self.populations2[state]*entropy2

    def calculate_all_jsd(self):
        """Calculates the Jennsen Shannon Divergence for all states"""
        jsd = [] 
        for state in range(self.sparse2.shape[0]):
            counts1, counts2 = self.prepare_count_arrays(state)
            jsd_for_state = self.calculate_jsd(state, counts1, counts2)
            jsd.append(jsd_for_state)
        return jsd

    def calculate_jsd_var(self, state,  n_bootstraps=50):
        """Calculate the variance of the JSD for "state" using 
           the bootstrap method with "n_boostraps" boostraps """
        counts1, counts2 = self.prepare_count_arrays(state)
        total_counts1 = np.sum(counts1)
        total_counts2 = np.sum(counts2)
 
        if total_counts1 == 0:      #If total counts is zero, 
            prob1 = counts1         #all probabilities are zero, which will be
        else:                       #equivalent to the count matrix.
            prob1 = np.divide(counts1, total_counts1)
        if total_counts2 == 0:              
            prob2 = counts2   
        else:
            prob2 = np.divide(counts2, total_counts2)

        jsds = []

        #Draw counts from a multinomial distribution
        sampled_counts1 = np.random.multinomial(total_counts1, prob1,  
                                                size = n_bootstraps)
        sampled_counts2 = np.random.multinomial(total_counts2, prob2,
                                                size = n_bootstraps) 
        for trial in range(n_bootstraps):
            jsd = self.calculate_jsd(state, sampled_counts1[trial,:],
                                     sampled_counts2[trial,:])
            jsds.append(jsd)
        avg_populations = ((self.populations1[state]+self.populations2[state])/2.0)**2
        return avg_populations*np.array(jsds).var()

    def calculate_all_jsd_var(self, nbootstraps=50):
        variances = []
        for istate in range(self.sparse2.shape[0]):
            variance = self.calculate_jsd_var(istate, nbootstraps) 
            if (variance < 10**-12):
                variances.append(0.0)
            else:
                variances.append(variance)
        return variances

if __name__ == '__main__':

    usage = "Surprisal.py OLD_SPARSE_MATRIX NEW_SPARSE_MATRIX"

    if len(sys.argv) < 3:
        print(usage)

    sparse1 = mmread(sys.argv[1])
    sparse2 = mmread(sys.argv[2])
    obj = SurprisalCalculator(sparse1, sparse2) 
    
    print "state\tsurprisal"

    for istate in range(sparse1.shape[0]):
        counts1, counts2 = obj.prepare_count_arrays(istate)
        surprisal = obj.calculate_surprisal(counts1, counts2, normalized=True)
        print "%d\t%e"%(istate, surprisal) 

