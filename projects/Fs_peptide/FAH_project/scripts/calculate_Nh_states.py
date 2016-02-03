import os,sys
import numpy as np


mapid = 40

traj_log_fn = '../../trajectories.log'
with open(traj_log_fn,'r') as logfn:
    trajlog = logfn.readlines()

all_Nhs = np.load("../../Number_helical_residues_all_data.npy")

for p_id in range(6383,6391):
    Nhs_state = []
    data_dir = "Data-{}-macro{}".format(p_id,mapid)
    assignment_fn = 'Assignments-{}.fixed.Map{}.npy'.format(p_id,mapid)
    a = np.load(assignment_fn)
    a =np.concatenate(a)
    project_index = [i for i,j in enumerate(trajlog) if str(p_id) in j]
    Nhs_project = all_Nhs[project_index]
    Nhs_project = np.concatenate(Nhs_project)
    for stateid in range(mapid):
        Nhs_state.append(np.mean(Nhs_project[np.where(a==stateid)]))
    output_fn = os.path.join(data_dir,"Nhs_state.dat")
    np.savetxt(output_fn,Nhs_state)
    print "Saved:{}".format(output_fn)

