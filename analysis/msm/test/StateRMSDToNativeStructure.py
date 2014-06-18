#----------------------------
#Written by Guangfeng Zhou, October,15,2012
#Dr.Voelz Lab
#Temple University
#----------------------------

import os,sys
sys.path.append('/Users/tud51931/scripts/gfzhou')
import HelixCoilTools as hct
import numpy as np
import copy
from scipy import savetxt
from msmbuilder import Serializer

cutoff = 3.0
metrics = 'rmsd'

if metrics.lower() == 'dihedral':
    Path = "/Users/tud51931/projects/MSM/msm/ff03-dihedralhybrid/"
    metrics = 'Dihedral'
elif metrics.lower() == 'rmsd':
    Path = "/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/"
    metrics = 'RMSD'
Path = os.path.join(Path,'%sCluster%0.1f'%(metrics,cutoff))   
    
AssignmentFile = os.path.join(Path,"Data","Assignments.h5")
A = Serializer.LoadFromHDF(AssignmentFile)
StateAssignment = hct.get_StatesAssignments(AssignmentFiles = A)
RMSDFile = os.path.join(Path,"Data","RMSD.h5")
RMSD = Serializer.LoadFromHDF(RMSDFile)
rmsd_allstates = {}
for state in StateAssignment.keys():
    rmsd_singlestate = []
    for trajid in StateAssignment[state].keys():
        rmsd_singlestate += list(RMSD['Data'][int(trajid)][StateAssignment[state][trajid]])
    rmsd_allstates[int(state)] = rmsd_singlestate

maxstatelength = max([len(i) for i in rmsd_allstates.values()])
StateRMSDs = copy.deepcopy(RMSD)
StateRMSDs['Data'] = -1*np.ones((len(rmsd_allstates),maxstatelength))
for state in rmsd_allstates.keys():
    statelength = len(rmsd_allstates[state])
    StateRMSDs['Data'][state][:statelength] = rmsd_allstates[state][:]
fn = "StateRMSDs_%sCluster%0.1f.h5"%(metrics,cutoff)
Serializer.SaveToHDF(StateRMSDs,fn)
print "Save to %s"%fn
print "Done."

    
