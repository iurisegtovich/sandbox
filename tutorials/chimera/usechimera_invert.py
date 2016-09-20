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
    inpdb = '/Users/tud51931/projects/QM/albocycline/conformations/%d.pdb'%i
    outpdb = inpdb.replace('conformations','conformations_12epi')
    rc("open %s"%inpdb)
    rc("invert #0@C6")
    rc("invert #0@C5")
    rc("write format pdb 0 %s"%outpdb)
    rc("close 0")
