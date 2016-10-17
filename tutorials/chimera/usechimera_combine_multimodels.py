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
    for j in range(5):
        recptor = '../../../../MurA-MSM-mol2/state%d-%d.aligned.mol2'%(i,j)
        ligand = '../../MurA-albo-dock/state%d-%d/rigid-xtal1_scored.mol2'%(i,j)
        outfile = 'state%d-%d-albo1.pdb'%(i,j)
        rc("open %s"%recptor)
        rc("open %s"%ligand)
        rc("combine #0,1.1")
        rc("write format pdb 2 %s"%outfile)
        rc("close all")
