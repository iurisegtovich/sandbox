import os,sys
import numpy as np
from scipy import loadtxt,savetxt
from msmbuilder.msm import MarkovStateModel

msm = MarkovStateModel(lag_time=50,n_timescales=10)
scores = []
#for pid in range(6383,6391):
#    for i in [100]:
#        assignment = np.load("Assignments-%d.fixed.Map%d.BACE+.npy"%(pid,i))
#        msm.fit(assignment)
#        scores.append([pid,i,msm.score_])
#output_fn = "GMRQ_scores_BACE+_MacroMSMs_extra_per_sequence.txt"
for i in [100]:
    assignment = np.load("Assignments.fixed.Map%d.BACE+.npy"%(i))
    msm.fit(assignment)
    scores.append([i,msm.score_])

output_fn = "GMRQ_scores_BACE+_MacroMSMs_extra.txt"

if os.path.exists(output_fn):
    print "%s exists! Exit!"%output_fn
else:
    savetxt(output_fn,scores,fmt="%d %f")

    
