import sys
import mdtraj
import numpy as np


#fn_pdb = '/Users/tud51931/projects/Helix-coil/data/Ala5/Ala5.pdb'
fn_pdb = sys.argv[1]

pdb = mdtraj.load_pdb(fn_pdb)
phi_indices = mdtraj.compute_phi(pdb)[0]
psi_indices = mdtraj.compute_psi(pdb)[0]

np.savetxt('Ala5_phi_indices.ndx',phi_indices,fmt='%d')
np.savetxt('Ala5_psi_indices.ndx',psi_indices,fmt='%d')

