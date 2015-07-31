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
            print "Matrix %d has a different shape!"%i
            raise




class SurprisalAnalysis:

    def __init__(self,*matrices,**kwargs):

        if len(matrices) < 2:
            raise("At least two count matrices!")
        check_matrices_shape(*matrices)

        self.tCounts_list = []
        for tCounts in matrices:
            self.tCounts_list.append(tCounts)

        self.n_matrices = len(args)



