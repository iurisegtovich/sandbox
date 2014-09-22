import os,sys
from msmbuilder import Project
import mdtraj as md
from mdtraj import io
import numpy as np

project = Project.load_from("ProjectInfo-RRR.yaml")
Rgs = -1*np.ones((project.n_trajs,max(project.traj_lengths)))

for i in range(project.n_trajs):
    t = project.load_traj(i)
    rg = md.compute_rg(t)
    Rgs[i][:len(rg)]=rg

io.saveh('Rgs-RRR.h5',Rgs)
