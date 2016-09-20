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
    outpdb = '%d.fixed.preminimized.pdb'%i
    rc("open %s"%inpdb)
    rc("minimize nogui #0@C11@H7@O3@C12@H2@H4@C19@H191@H192@H193")
    rc("write format pdb 0 %s"%outpdb)
    rc("close 0")
