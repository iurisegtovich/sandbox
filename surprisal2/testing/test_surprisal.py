import os,sys
from scipy.io import mmread
sys.path.append("../../src/")
from Surprisal_GZhou import SurprisalAnalysis

matrix1 = mmread("../wt-bb-cb-2000-2ev-t-20macro-tCounts.mtx").tolil()
matrix2 = mmread("../tz4-bb-cb-2000-2ev-t-20macro-tCounts.mtx").tolil()

obj = SurprisalAnalysis(matrix_type = "sparse",normalize = "counts",var_method="analytical")
obj.calculate_surprisals(matrix1,matrix2)
obj.calculate_surprisals_var(matrix1,matrix2)
obj.calculate_JSDs(matrix1,matrix2)
print obj.surprisals_
print obj.surprisals_var_
print obj.surprisal_weights_
print obj.JensenShannonDivergence_
