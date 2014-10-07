from mdtraj import io
import argparse
import numpy as np

def mean_Nv_state(assignment,Nv):
    nv_state = np.zeros(assignment.max()+1)
    for stateid in range(assignment.max()+1):
        nv_state[stateid] = Nv[assignment==stateid].mean()
    return nv_state



assign_file = "/Volumes/Guangfeng/Fs-peptide/Fs-ff03-owlsnest/HelixCoil/Data/Assignments.h5"
Nv_file = "Nv.h5"

assignment = io.loadh(assign_file,'arr_0')
Nv = io.loadh(Nv_file,'arr_0')

Mean_Nv_state = mean_Nv_state(assignment,Nv)
print Mean_Nv_state
np.savetxt('mean_Nv_state.dat',Mean_Nv_state)
