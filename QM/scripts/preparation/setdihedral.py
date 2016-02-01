import __main__
__main__.pymol_argv = ['pymol','-qc']
import pymol
import numpy as np
from pymol import cmd
pymol.finish_launching()

def _set_dihedral(index,angle):
    for i in range(4):
        cmd.select('a%d'%i, 'id %d'%index[i])
    cmd.set_dihedral('a0','a1','a2','a3',angle)

def read_indices(filename):
    f = open(filename,'r')
    lines = f.readlines()
    f.close()
    indices = []
    for line in lines:
        indices.append(line.replace('\n','').split(','))
    return indices

def Set_Dihedral(indices,angle):
    for i in indices.shape[0]:
        _set_dihedral(indices[i],angle)
