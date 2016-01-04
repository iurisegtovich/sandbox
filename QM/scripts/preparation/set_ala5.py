import os,sys
import numpy as np
import setdihedral

phi_ndx = np.loadtxt('/Users/tud51931/projects/Helix-coil/data/Ala5/prepare/Ala5_phi_indices.ndx')
psi_ndx = np.loadtxt('/Users/tud51931/projects/Helix-coil/data/Ala5/prepare/Ala5_psi_indices.ndx')
input_pdb = '/Users/tud51931/projects/Helix-coil/data/Ala5/prepare/Ala5.pdb'
cmd.load(input_pdb,'mol')

#pymol atom id start from 1
phi_ndx = phi_ndx + 1
psi_ndx = psi_ndx + 1

#alpha-helix
Phi_h = -60.0
Psi_h = -45.0
#coil
Phi_c = 180
Psi_c = 180

for i in range(32):
    for j in range(5):
        try:
            if bin(i)[-j-1]=='1':
                setdihedral._set_dihedral(phi_ndx[j],Phi_h)
                setdihedral._set_dihedral(psi_ndx[j],Psi_h)
            elif bin(i)[-j-1]=='0':
                setdihedral._set_dihedral(phi_ndx[j],Phi_c)
                setdihedral._set_dihedral(psi_ndx[j],Psi_c)
        except IndexError:
            setdihedral._set_dihedral(phi_ndx[j],Phi_c)
            setdihedral._set_dihedral(psi_ndx[j],Psi_c)

    fn_outpdb = 'Ala5_%d.pdb'%i
    cmd.save(fn_outpdb,'mol')
