#Assign the snapshot according to helix coil sequence.
import os,sys
from msmbuilder import Project
from mdtraj import io
import mdtraj
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-p','--Project',help='path to the project.yaml file',required=True)
parser.add_argument('-o','--Output',help='Assignment output file',default='Assignments.h5')
args = parser.parse_args()



def is_helical_peptide(phiangle, psiangle):

        phirange = [-60.0-30.0, -60.0+30.0]
        psirange = [-47.0-30.0, -47.0+30.0]

        if phiangle > phirange[0] and phiangle < phirange[1]:
            if psiangle > psirange[0] and psiangle < psirange[1]:
                return True
        else:
            return False


def ConvertDihedralsToArray(phi,psi):
    HCarray = np.zeros((phi.shape))
    for i in range(len(phi)):
        for j in range(len(phi[i])):
            if is_helical_peptide(phi[i,j],psi[i,j]):
                HCarray[i][j] = 1
            else:
                HCarray[i][j] = 0
    return HCarray

def count_n_helices(HCarray):

    return HCarray.sum(1)





project = Project.load_from(args.Project)
assignments = -1*np.ones((project.n_trajs,max(project.traj_lengths)))

for trajid in range(project.n_trajs):
    print "Working on: %s"%project.traj_filename(trajid)
    traj = project.load_traj(trajid)
    phi = mdtraj.compute_phi(traj)[1]*360/(2*np.pi)
    psi = mdtraj.compute_psi(traj)[1]*360/(2*np.pi)
    HCarray = ConvertDihedralsToArray(phi,psi)
    assignments[trajid][:traj.n_frames] = count_n_helices(HCarray)



io.saveh(args.Output,assignments)