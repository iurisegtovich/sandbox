import __main__
__main__.pymol_argv = ['pymol','-qc']
import pymol
import numpy as np
from pymol import cmd
pymol.finish_launching()

def _set_dihedral(atoms,angle,mode = "index"):
    if mode == "name":
        for i in range(4):
            cmd.select('a%d'%i, 'name %s'%atomnames[i])
    elif mode == "index":
        for i in range(4):
            cmd.select('a%d'%i, 'index %d'%atoms[i])
    cmd.set_dihedral('a0','a1','a2','a3',angle)

def read_indices(filename):
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    indices = []
    for line in lines:
        indices.append(line.replace('\n','').split(','))
    return indices

def Set_Dihedral(indices,angle,mode="index"):
    for index in indices:
        if isinstance(index,list):
            _set_dihedral(index,angle,mode)
        else:
            _set_dihedral(indices,angle,mode)

phi_ndx = [12,1,3,14]
psi_ndx = [1,3,14,2]
omega_ndx = [13,12,1,3]
chi1_ndx = [5,4,1,3]
chi2_ndx = [7,6,4,5]

input_pdb = '../template/P5-cis.pdb'

bb_conf = "cis-a"

if bb_conf == "cis-a":
    omega_angle = 0.0
    phi_angle = -90.0
    psi_angle = 180.0
elif bb_conf == "trans-a":
    omega_angle = 180.0
    phi_angle = -90.0
    psi_angle = 180.0
elif bb_conf == "trans-b":
    omega_angle = 180.0
    phi_angle = -130.0
    psi_angle = 80.0

cmd.load(input_pdb,'mol')
i = 0
Set_Dihedral(phi_ndx,phi_angle)
Set_Dihedral(psi_ndx,psi_angle)
Set_Dihedral(omega_ndx,omega_angle)
angles = []
for chi1 in np.arange(0,360,30.):
    for chi2 in np.arange(0,180,30.):
        Set_Dihedral(chi1_ndx,chi1)
        Set_Dihedral(chi2_ndx,chi2)
        output_pdb = 'chi1chi2-%s-%d.pdb'%(bb_conf,i)
        cmd.save(output_pdb,'mol')
        i += 1
        print "Wrote: %s"%output_pdb
        angles.append([chi1,chi2])

angle_fn = 'chi1chi2_%s_angles_index.dat'%bb_conf
np.savetxt(angle_fn,angles)

"""
cmd.load(input_pdb,'mol')
Set_Dihedral(PhiIndices,PhiAngle)
Set_Dihedral(PsiIndices,PsiAngle)
Set_Dihedral(OmegaIndices,OmegaAngle)
Set_Dihedral(ChiIndices,ChiAngle)
cmd.save(output_pdb,'mol')
"""
