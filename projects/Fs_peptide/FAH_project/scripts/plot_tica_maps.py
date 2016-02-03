import os,sys
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm
from matplotlib.colors import LogNorm
from sklearn.externals import joblib
import mdtraj as md

tica_data = np.load('../../ticas_n_8.npy')
tica_data = np.concatenate(tica_data)

plt.figure()
plt.hist2d(tica_data[:,0],tica_data[:,1],bins=500,norm=LogNorm())
plt.xlabel("First tIC")
plt.ylabel("Second tIC")
plt.title("tICA Heatmap(log color scale)")
#plt.savefig("tica_all.eps")
#plt.show()

msm_fn = 'MSMs-combined-macro40/MSMs-combined.pkl'
traj_log_fn = '../../trajectories.log'
assignment_fn = 'Assignments.fixed.Map40.npy'

msm = joblib.load(msm_fn)
a = np.load(assignment_fn)
selected_pairs_by_state = msm.draw_samples(a,50)
tica_data = np.load('../../ticas_n_8.npy')

color=iter(cm.rainbow(np.linspace(0,1,40)))

for state_id in range(selected_pairs_by_state.shape[0]):
    c = next(color)
    for i,(traj_id,frame_id) in enumerate(selected_pairs_by_state[state_id]):
        plt.plot(tica_data[traj_id][frame_id][0],tica_data[traj_id][frame_id][1],marker='o',c=c)
plt.show()

#PDBs_dir = 'PDBs-combined'
#if not os.path.exists(PDBs_dir):
#    os.makedirs(PDBs_dir)

#with open(traj_log_fn,'r') as logfn:
#    trajlog = logfn.readlines()
#for state_id in range(selected_pairs_by_state.shape[0]):
#    for i,(traj_id,frame_id) in enumerate(selected_pairs_by_state[state_id]):
#        p_id = trajlog[traj_id][:-1].split('/')[-2][-4:]
#        pdb_fn = '../../Gen0-%s.pdb'%p_id
#        pdb = md.load_pdb(pdb_fn)
#        traj = md.load(trajlog[traj_id][:-1],top=pdb)
#        output_fn = os.path.join(PDBs_dir,"State%d-%d.pdb"%(state_id,i))
#        traj[frame_id].save_pdb(output_fn)

