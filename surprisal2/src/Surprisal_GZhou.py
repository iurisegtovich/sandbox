#/usr/bin/env python

import sys
import os 
import numpy as np
from scipy.io import mmread, mmwrite
from scipy.sparse import isspmatrix_coo	

#from sklearn.preprocessing import normalize
#from msmbuilder.MSMLib import estimate_transition_matrix, mle_reversible_count_matrix
#from msmbuilder import msm_analysis     # TODO: check version?
#from Reversible_MLE_with_Guess import *

def check_matrices_shape(*matrices):

    s = matrices[0].shape
    for i in range(len(matrices)):
        if matrices[i].shape != s:
            raise TypeError("Matrix %d has a different shape!"%i)
    return s

class SurprisalAnalysis:

    def __init__(self, **kwargs):
        """
        Calculates surprisal for arbitrary sets of counts. If normalized, the
        surprisal value is divided by the total number of counts from that state.

        INPUTS
        c1, c2, ....   any number of numpy counts

        OPTIONS
        normalize   If True, will divide the result by the total numnber of counts
        """


        self.normalize = "counts"
        for k,v in kwargs.iteritems():
            if k == 'normalize':
                self.normalize = v

    def calculate_surprisal_state(self, state_id, *args):

        """
        Calculates surprisal for arbitrary sets of counts. If normalized, the
        surprisal value is divided by the total number of counts from that state.
        INPUTS
        c1, c2, ....   any number of numpy counts
        OPTIONS
        normalize   If True, will divide the result by the total numnber of counts
        """

        if len(args) < 2:
            raise ValueError(' must be given at least two count arrays!')

        # compute combined counts
        c_comb = args[0][state_id] + args[1][state_id]
        if len(args) > 2:
            for i in range(2, len(args)):
                c_comb += args[i][state_id]


        # compute count totals
        totals = [c[state_id].sum() for c in args]
        total_comb = c_comb.sum()
        if total_comb != 0:
            surprisal = float(total_comb)

        if self.normalize == "counts":

            # Voelz, V. A., Elman, B., Razavi, A. M., & Zhou, G. (2014).
            # Surprisal Metrics for Quantifying Perturbed Conformational Dynamics in Markov State Models.
            # Equation (7)
            surprisal = H(c_comb)
            for i in range(len(args)):
                surprisal -= totals[i]/total_comb*H(args[i][state_id])

        elif self.normalize == "mle":
            # Voelz, V. A., Elman, B., Razavi, A. M., & Zhou, G. (2014).
            # Surprisal Metrics for Quantifying Perturbed Conformational Dynamics in Markov State Models.
            # Equation (14) in the parentheses.
            print "Not implemented yet"

        return surprisal


    def calculate_surprisal(self, *matrices):

        self.surprisal_total_ =[]
        matrix_shape = check_matrices_shape(*matrices)

        for state_id in range(matrix_shape[0]):
            self.surprisal_total_.append(self.calculate_surprisal_state(state_id,*matrices))



def counts2probs(c):
    # convert to float if necessary
    if c.dtype.kind == 'i':
        c = c.astype('float')
    return c/c.sum()


def H(p, normalize=True):
    """Returns the entropy H = \sum_i - p_i ln p_i of a distribution.
    INPUTS
    p	        numpy array of values (discrete counts work too)
    OPTIONS
    normalize 	If True, will normalize the input values so that \sum_i p_i = 1
    """

    if normalize:
        p = counts2probs(p)

    # non-zero entries only
    Ind = (p>0)

    return np.dot(-p[Ind], np.log(p[Ind]))

def H_cross(p,q,normalize=True):

    if normalize:
        p = counts2probs(p)
        q = counts2probs(q)

    # non-zero entries only
    Ind = (p>0)
    return np.dot(-p[Ind], np.log(q[Ind]))

def H_var(c):
    """Estimate the variance of the entropy H due to finite sampling.
    INPUT
    c     	a row of transition counts
    OUTPUT
    var_H	an estimate of the variance of H(c)
    """

    # build MVN/multinomial covariance matrices for c1 and c2
    V = cov_multinomial_counts(c)

    # compute the vector of sensitivities q = [dH/dnj ... ]
    n = c.sum()
    print 'c =', c, 'n =', n
    q = H(c)*np.ones( c.shape )
    for i in range(c.shape[0]):
        if c[i] > 0:
            q[i] -= np.log( float(c[i])/n )
    q = q/n

    # return q^T V q
    return np.dot(q, V.dot(q))


