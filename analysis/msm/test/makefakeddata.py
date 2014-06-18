#makefakeddata.py
#Make faked data to test the Surprisal value calculation

import os,sys
import numpy as np
from scipy.io import mmread
from msmbuilder import Serializer


def draw_index(probs, n_picks=1, UseFastMethod=True):
    """Draw a number (or many numbers, controlled by n_picks), weighted by the probabilities probs."""
    if UseFastMethod:
        t = np.cumsum(probs)
        s = sum(probs)
        return np.searchsorted(t,np.random.rand(n_picks)*s)

tcounts = mmread('/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster4.2/lagtime50/tCounts.UnMapped.mtx')
Assignment = Serializer.LoadFromHDF('/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster4.2/lagtime50/Assignments.Fixed.h5')
trajnum = 100
frames = Assignment['Data'].shape[1]
a = Serializer()
a['Data'] = -1*np.ones((trajnum,frames))
for traj in range(trajnum):
    print '%d of %d Trajectories'%(traj,trajnum)
    startstate = 126
    a['Data'][traj,0] = startstate
    for step in range(1,frames):
        probs = tcounts.data[tcounts.row == startstate]/sum(tcounts.data[tcounts.row == startstate])    
        a['Data'][traj,step] = tcounts.col[tcounts.row == startstate][draw_index(probs)[0]]
        startstate = a['Data'][traj,step]

a.SaveToHDF('Assignments.rmsd4.2.faked.h5')



    