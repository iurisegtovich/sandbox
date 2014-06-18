import os,sys
import chimera
import BuildStructure
from chimera import runCommand as rc


for i in range(100):
    inpdb = '../qm_minimized_albocycline_pdbs/%d.pdb'%i
    opened = chimera.openModels.open(inpdb)
    mol = opened[0]
    #Atoms[21] is H5. Changed H5 to CH3, 5-methyl-albocycline .
    BuildStructure.changeAtom(mol.atoms[21],'C',4,4)
    outpdb = inpdb.replace('qm_minimized_albocycline','5-methyl')
    rc("write format pdb 0 %s"%outpdb)
    rc("close 0")
