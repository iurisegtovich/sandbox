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

ref_pdb = ""
rc("open %s"%ref_pdb)
for i in range(100):
    inpdb = '%d.fixed.pdb'%i
    outpdb = '%d.fixed.preminimized.pdb'%i
    rc("open %s"%inpdb)
    rc("matchmaker #0 #1")
    rc("write format pdb 1 %s"%outpdb)
    rc("close 1")
