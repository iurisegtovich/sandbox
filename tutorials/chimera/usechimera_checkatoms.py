import os,sys
import chimera
import BuildStructure
from chimera import runCommand as rc

inpdb = '../qm_minimized_albocycline_pdbs/0.pdb'
opened = chimera.openModels.open(inpdb)
mol = opened[0]
for i in range(len(mol.atoms)):
    print i, mol.atoms[i].name
