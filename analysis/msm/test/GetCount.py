import sys, os
import scipy.io
from msmbuilder import MSMLib
from msmbuilder.io import loadh,saveh

try:
    Assignments=loadh("%s"%(sys.argv[1]),'arr_0').astype(int)
except KeyError:
    Assignments=loadh("%s"%(sys.argv[1]),'Data').astype(int)
NumStates = max(Assignments.flatten()) + 1
LagTime = int(sys.argv[2])
Counts = MSMLib.get_count_matrix_from_assignments(Assignments, n_states=NumStates, lag_time=LagTime, sliding_window=True)
scipy.io.mmwrite('%s'%(sys.argv[3]), Counts)
