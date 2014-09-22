from mdtraj import io
from collections import Counter
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-a','--Assignments',help='Path to the Assignments.h5 file',required=True)
parser.add_argument('-o','--Output',help='Stationary population output file name',default='Populations.dat')
args=parser.parse_args()

a = io.loadh(args.Assignments,'arr_0')
a = a.reshape(1,-1)
c = Counter(a[0])
population = np.zeros(np.max(c.keys())+1)
for stateid in range(np.max(c.keys())+1):
    population[stateid] = c[stateid]

population = population/population.sum()
print "Stationary Populaiton",population
np.savetxt(args.Output,population)



