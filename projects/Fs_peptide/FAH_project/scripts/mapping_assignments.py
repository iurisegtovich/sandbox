import os,sys
import numpy as np
from scipy import loadtxt

def map_assignments(mapping,assignments):
    fixed_assignments = []
    for i in range(assignments.shape[0]):
        for j in range(len(assignments[i])):
            assignments[i][j] = mapping[assignments[i][j]]
        fixed_assignments.append(assignments[i])
    return fixed_assignments


def case1():
    ids = [30,40,50]

    for map_id in ids:
        mapping = loadtxt("Output_BACE/map%d.dat"%map_id)
        assignments = np.load("Assignments.npy")
        fixed_assignments = map_assignments(mapping,assignments)
        fixed_assignments_fn = "Assignments.fixed.Map%d.npy"%(map_id)
        np.save(fixed_assignments_fn,fixed_assignments)
        print "Wrote: %s"%fixed_assignments_fn

def case2():
    map_id = 40
    mapping = loadtxt("Output_BACE/map%d.dat"%map_id)
    for p_id in range(6383,6391):
        assignments = np.load("Assignments-%d.npy"%p_id)
        fixed_assignments = map_assignments(mapping,assignments)
        output_fn = "Assignments-%d.fixed.Map%d.npy"%(p_id,map_id)
        np.save(output_fn,fixed_assignments)
        print "Wrote: %s"%output_fn
def case3():
    map_id = 100
    mapping = loadtxt("/Users/tud51931/git/voelzlab/baceplus/src/OUT_BACE+/map%d.dat"%map_id)
    for p_id in range(6383,6391):
        assignments = np.load("Assignments-%d.npy"%p_id)
        fixed_assignments = map_assignments(mapping,assignments)
        output_fn = "Assignments-%d.fixed.Map%d.BACE+.npy"%(p_id,map_id)
        np.save(output_fn,fixed_assignments)
        print "Wrote: %s"%output_fn
    assignments = np.load("Assignments.npy")
    fixed_assignments = map_assignments(mapping,assignments)
    fixed_assignments_fn = "Assignments.fixed.Map%d.BACE+.npy"%(map_id)
    np.save(fixed_assignments_fn,fixed_assignments)
    print "Wrote: %s"%fixed_assignments_fn



if __name__ == '__main__':
    case3()

        

    
    
