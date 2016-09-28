#------------------------------------------
#Written by Guangfeng Zhou
#Dr.Voelz Lab
#Chemistry Department
#Temple University
#October,25,2012
#-------------------------------------------
import os,sys
import chimera
from chimera import runCommand as rc

for i in range(100):
    inpdb = '%d.fixed.pdb'%i
    outmol2 = '%d.fixed.mol2'%i
    rc("open %s"%inpdb)
    rc("del #0&H")
    rc("addh")
    rc("addcharge")
    rc("write format mol2 0 %s"%outmol2)
    rc("close 0")
