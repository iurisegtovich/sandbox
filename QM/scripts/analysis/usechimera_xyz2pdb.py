#-------------------------------------------
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
    xyz = 'outputs/%d.xyz'%i
    pdb = xyz.replace('.xyz','.pdb')
    rc("open %s"%xyz)
    rc("write format pdb 0 %s"%pdb)
    rc("close 0")
