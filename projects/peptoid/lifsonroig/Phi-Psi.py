import os,sys
import mdtraj as md
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

dir_path = "/home/tue91994/PROJECTS/Peptoids/SPE_Nmers/SPE_7mer"
pdb_fn = os.path.join(dir_path,"SPE7.pdb")
pdb = md.load_pdb(pdb_fn)
trajs = []

for i in range(1):
    traj_fn = os.path.join(dir_path,"traj%d.xtc"%i)
    trajs.append(md.load(traj_fn,top=pdb))
phi = md.compute_phi(trajs)
psi = md.compute_psi(trajs)

plt.figure()
plt.hist2d(phi,psi,bins=100,norm=LogNorm())
plt.colorbar()
plt.savefig("phi-psi.pdf")

