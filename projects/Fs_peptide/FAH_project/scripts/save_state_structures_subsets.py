import os,sys
import numpy as np
import mdtraj as md
from sklearn.externals import joblib

traj_log_fn = '../../trajectories.log'
with open(traj_log_fn,'r') as logfn:
    trajlog = logfn.readlines()
for p_id in range(6383,6391):
    msm_fn = 'MSMs-%d-macro40/MSMs-%d-macro40.pkl'%(p_id,p_id)
    assignment_fn = 'Assignments-%d.fixed.Map40.npy'%p_id
    trajlog_id = [j for i,j in enumerate(trajlog) if str(p_id) in j]

    pdb_fn = '../../Gen0-%s.pdb'%p_id
    pdb = md.load_pdb(pdb_fn)
    msm = joblib.load(msm_fn)
    a = np.load(assignment_fn)
    selected_pairs_by_state = msm.draw_samples(a,5)

    PDBs_dir = 'PDBs-%d'%p_id
    if not os.path.exists(PDBs_dir):
        os.makedirs(PDBs_dir)

    for state_id in range(selected_pairs_by_state.shape[0]):
        for i,(traj_id,frame_id) in enumerate(selected_pairs_by_state[state_id]):
            traj = md.load(trajlog_id[traj_id][:-1],top=pdb)
            output_fn = os.path.join(PDBs_dir,"State%d-%d.pdb"%(state_id,i))
            traj[frame_id].save_pdb(output_fn)
