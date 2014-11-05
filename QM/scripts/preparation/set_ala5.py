import os,sys
import numpy as np
import usepymol_setdihedral as setdihedral

phi_ndx = np.loadtxt('')
psi_ndx = np.loadtxt('')
input_pdb = ''
cmd.load(input_pdb,'mol')

#alpha-helix
Phi_h = -60.0
Psi_h = -45.0
#coil
Phi_c = 10.0
Psi_c = 10.0

for i in range(32):
    for j in range(5):
        if bin(i)[-j]=='1':
            setdihedral.Set_Dihedral(phi_ndx[j],Phi_h)
            setdihedral.Set_Dihedral(psi_ndx[j],Psi_h)
        elif bin(i)[-j]=='0':
            setdihedral.Set_Dihedral(phi_ndx[j],Phi_c)
            setdihedral.Set_Dihedral(psi_ndx[j],Psi_c)

    fn_outpdb = 'Ala5_%d.pdb'%i
    cmd.save(fn_outpdb,'mol')
