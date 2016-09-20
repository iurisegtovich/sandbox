import os,sys
import chimera
import BuildStructure
from chimera import runCommand as rc

def findindex(atomname):
    for i in range(len(mol.atoms)):
        if mol.atoms[i].name.lower() == atomname.lower():
            return i

for i in range(100):
    inpdb = '../qm_minimized_albocycline_pdbs/%d.pdb'%i
    opened = chimera.openModels.open(inpdb)
    mol = opened[0]
    index = findindex('H3')
    #rc('delete #0@H8@H9@H10')
    BuildStructure.changeAtom(mol.atoms[index],'C',4,4)
    outpdb = '%d.pdb'%i
    rc("write format pdb 0 %s"%outpdb)
    rc("close 0")
