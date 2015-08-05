import os,sys
from scipy.io import mmread
sys.path.append("../src/")
from Surprisal import SurprisalAnalysis
import numpy as np

matrices = []
for pid in range(6383,6391):
    matrix = mmread("Fs-%d-tCounts-macro40.mtx"%pid).tolil()
    matrices.append(matrix)

obj = SurprisalAnalysis(matrix_type = "sparse",normalize = "counts",var_method="analytical")
obj.calculate_surprisals(*matrices)
obj.calculate_surprisals_var(*matrices)
obj.calculate_JSDs(*matrices)
print obj.surprisals_
print obj.surprisals_var_
print obj.surprisal_weights_
print obj.JensenShannonDivergence_
print np.argsort(obj.JensenShannonDivergence_)
