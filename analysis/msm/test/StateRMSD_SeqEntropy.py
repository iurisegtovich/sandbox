#----------------------------
#October,17,2012
#Guangfeng Zhou
#Dr.Voelz Lab
#Room 100/102, Beury Hall
#Temple University

import sys,os
sys.path.append('/Users/tud51931/scripts/gfzhou')
from msmhcanalysis import SequenceEntropy_states
from msmbuilder import Serializer
import matplotlib.pyplot as plt
import numpy as np

stateentropy = SequenceEntropy_states('HCstrings_states_Dihedral5.2.txt')
rmsdfile = Serializer.LoadFromHDF('StateRMSDs_DihedralCluster5.2.h5')
RMSD = np.ma.array(rmsdfile['Data'],mask=[rmsdfile['Data']==-1])
statermsd = RMSD.mean(1)

plt.figure()
plt.plot(statermsd,stateentropy,'.')
plt.title('StateSequecneEntropy versus StateRMSD')
plt.ylabel('StateSequenceEntropy')
plt.xlabel('StateRMSD(nm)')
plt.savefig('seqentropy_statermsd_dihedralcluster.png')
#plt.show()


