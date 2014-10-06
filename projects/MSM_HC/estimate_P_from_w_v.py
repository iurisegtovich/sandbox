import os,sys
sys.path.append('/Users/tud51931/voelzlab/analysis/LifsonRoig/src')
from mdtraj import io
import mdtraj
from msmbuilder import Project
import LifsonRoigTools as LRTools
import numpy as np
from collections import Counter

v = 0.351
w = 1.615

assignment = io.loadh('/Volumes/Guangfeng/Fs-peptide/Fs-ff03-owlsnest/ff03-RRR-full-msm/Assignments.h5.RRR','arr_0')
project = Project.load_from('/Volumes/Guangfeng/Fs-peptide/Fs-ff03-owlsnest/ff03-RRR-full-msm/ProjectInfo.yaml')
c = Counter(assignment.reshape(1,-1)[0])
populations = np.zeros(np.max(c.keys())+1)

for i in range(project.n_trajs):
    print "Working on: %s"%project.traj_filename(i)
    traj = project.load_traj(i)
    hc = LRTools.LifsonRoigTools()
    hc.prepare(traj)
    hc.ConvertDihedralsToHCStrings(sequencetype='peptide',helixtype='alpha')
    hc.count_Helix_per_frame()
    for j in range(assignment[i].shape[0]):
        if assignment[i][j] != -1:
            populations[assignment[i][j]] = populations[assignment[i][j]] + (w**hc.helix_w_numbers_frame[j])*(v**hc.helix_v_numbers_frame[j])

# weigh populations by the counts
for i in range(populations.shape[0]):
    populations[i] = populations[i]/c[i]
populations = populations/populations.sum()
print populations
np.savetxt('Populations.LifsonRoig.dat',populations)

