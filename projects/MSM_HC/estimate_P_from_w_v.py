import os,sys
sys.path.append('/Users/tud51931/voelzlab/analysis/LifsonRoig/src')
from mdtraj import io
import mdtraj
from msmbuilder import Project
import LifsonRoigTools as LRTools
import numpy as np
from collections import Counter

l = np.loadtxt('/Users/tud51931/voelzlab/analysis/LifsonRoig/scripts/test_Fs_RRR_ff03Loglikelihood.dat')
w = np.loadtxt('/Users/tud51931/voelzlab/analysis/LifsonRoig/scripts/test_Fs_RRR_ff03w_params.dat')
v = np.loadtxt('/Users/tud51931/voelzlab/analysis/LifsonRoig/scripts/test_Fs_RRR_ff03v_params.dat')

I = np.argsort(l)
w_max = w[I[-1]]
v_max = v[I[-1]]

assignment = io.loadh('/Volumes/Guangfeng/Fs-peptide/Fs-ff03-owlsnest/ff03-RRR-full-msm/Assignments.h5.RRR','arr_0')
project = Project.load_from('/Volumes/Guangfeng/Fs-peptide/Fs-ff03-owlsnest/ff03-RRR-full-msm/ProjectInfo.yaml')
c = Counter(assignment.reshape(1,-1)[0])
populations = np.zeros(np.max(c.keys())+1)

def calculate_weight_frame(w_array,v_array,w_param,v_param):
    weight = 1.0
    for i in w_array*w_param:
        if i != 0:
            weight = weight*i
    for i in v_array*v_param:
        if i != 0:
            weight = weight*i
    return weight

for i in range(project.n_trajs):
    print "Working on: %s"%project.traj_filename(i)
    traj = project.load_traj(i)
    hc = LRTools.LifsonRoigTools()
    hc.prepare(traj)
    hc.ConvertDihedralsToHCStrings(sequencetype='peptide',helixtype='alpha')
    hc.calculate_w_v_array()

    for j in range(assignment[i].shape[0]):
        if assignment[i][j] != -1:
            weight = calculate_weight_frame(hc.w_array[j],hc.v_array[j],w_max,v_max)
            populations[assignment[i][j]] = populations[assignment[i][j]] + weight

# weigh populations by the counts
for i in range(populations.shape[0]):
    populations[i] = populations[i]/c[i]
populations = populations/populations.sum()
print populations
np.savetxt('Populations.LifsonRoig.new.dat',populations)

