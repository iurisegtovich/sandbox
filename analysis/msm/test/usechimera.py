#-------------------------------------------
#Written by Guangfeng Zhou
#Dr.Voelz Lab
#Chemistry Department
#Temple University
#October,25,2012
#-------------------------------------------
import chimera
from chimera import runCommand as rc

for i in range(0,1):
    rc("open gen%d.pdb"%i)
    rc("swapaa glu #%d:10"%i)
    rc("write format pdb %d gen%d.ERR.pdb"%(i,i))
rc("close all")
