import __main__
__main__.pymol_argv = ['pymol','-qc']
import pymol
import numpy as np
from pymol import cmd
pymol.finish_launching()

def _set_dihedral(atomnames,angle):
    for i in range(4):
        cmd.select('a%d'%i, 'name %s'%atomnames[i])
    cmd.set_dihedral('a0','a1','a2','a3',angle)

def read_indices(filename):
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    indices = []
    for line in lines:
        indices.append(line.replace('\n','').split(','))
    return indices

def Set_Dihedral(indices,angle):
    for index in indices:
        if isinstance(index,list):
            _set_dihedral(index,angle)
        else:
            _set_dihedral(indices,angle)

PhiIndices = read_indices('PhiIndices.dat')
PsiIndices = read_indices('PsiIndices.dat')
OmegaIndices = read_indices('OmegaIndices.dat')
ChiIndices = read_indices('ChiIndices-ACE-NDM.dat')

# Get the 3rd,4th residue indices
if 1:
    PhiIndices = PhiIndices[2:4]
    PsiIndices = PsiIndices[2:4]
    OmegaIndices = OmegaIndices[1:3]
    ChiIndices = ChiIndices[2:4]

input_pdb = '/Users/tud51931/projects/peptoid/Nspe5/Nspe5_ACE_NMD_144/starting_conf_minimization/output.pdb'
#input_pdb = './Nspe5_pymol.pdb'
#output_pdb = './Nspe5_pymol.pdb'
PhiAngle = -75.0
PsiAngle = 180.0
OmegaAngle = 0.0
ChiAngle = 90.0

#Cis/Trans AlphaD,C7beta angles. Order: Omega,Phi,Psi.
angles = [[0.0,75.0,180.0],[0.0,-75.0,180.0],[180.0,75.0,180.0],[180.0,-75.0,180.0],[180.0,130.0,-80.0],[180.0,-130.0,80.0]]
chi_angles = [-90.0,60.0]
cmd.load(input_pdb,'mol')
conf_arr = np.nan*np.ones((144,8))
i=0
for omega1,phi1,psi1 in angles:
    for chi1 in chi_angles:
        Set_Dihedral(PhiIndices[0],phi1)
        Set_Dihedral(PsiIndices[0],psi1)
        Set_Dihedral(OmegaIndices[0],omega1)
        Set_Dihedral(ChiIndices[0],chi1)
        for omega2,phi2,psi2 in angles:
            for chi2 in chi_angles:
                Set_Dihedral(PhiIndices[1],phi2)
                Set_Dihedral(PsiIndices[1],psi2)
                Set_Dihedral(OmegaIndices[1],omega2)
                Set_Dihedral(ChiIndices[1],chi2)
                output_pdb = '/Users/tud51931/projects/peptoid/Nspe5/Nspe5_ACE_NMD_144/144confs-ACE-NDM/Nspe5-ACE-NDM-%d.pdb'%i
                cmd.save(output_pdb,'mol')
                conf_arr[i][:]=omega1,phi1,psi1,chi1,omega2,phi2,psi2,chi2
                i += 1
                print "Wrote: %s"%output_pdb
np.savetxt('conf_arr.dat',conf_arr)
np.savetxt('/Users/tud51931/projects/peptoid/Nspe5/Nspe5_ACE_NMD_144/144confs-ACE-NDM/conf_arr.dat',conf_arr)

    
"""
cmd.load(input_pdb,'mol')
Set_Dihedral(PhiIndices,PhiAngle)
Set_Dihedral(PsiIndices,PsiAngle)
Set_Dihedral(OmegaIndices,OmegaAngle)
Set_Dihedral(ChiIndices,ChiAngle)
cmd.save(output_pdb,'mol')
"""
