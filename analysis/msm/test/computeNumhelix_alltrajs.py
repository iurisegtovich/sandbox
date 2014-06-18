import os,sys
import numpy as np
from msmbuilder import Trajectory
sys.path.append('~/scripts/gfzhou/')
import HelixCoilTools as hct
from scipy import savetxt
"""
This script is to get the number of helix from trajectories.
"""
datafile = "./numhelix_alltraj.txt"
if os.path.exists(datafile):
    print "%s already exists!"%datafile
    print "quit."
    sys.exit()

path ="/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/sourcedata/Trajectories"
numhelix_alltraj = -1*np.ones((100,8000),dtype=int)
for i in range(100):
    Trajfile = "%s/trj%d.lh5"%(path,i)
    if os.path.exists(Trajfile):
        T = Trajectory.LoadFromLHDF(Trajfile)
        print "Compute number of helix for %s"%Trajfile
        numhelix = hct.compute_numhelix_trajectory(T)
        numhelix_alltraj[i][:len(numhelix)] = numhelix[:]

print "Save data to %s"%datafile
savetxt(datafile,numhelix_alltraj)
print "Done." 
