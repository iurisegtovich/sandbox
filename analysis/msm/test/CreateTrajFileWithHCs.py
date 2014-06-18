import os,sys
sys.path.append('~/scripts/gfzhou/')
import HelixCoilTools as hct
from msmbuilder import Trajectory
"""
This script shows how to create new trj files with hc strings.
"""
path ="/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/sourcedata/Trajectories"
for i in range(0,100):
    Trajfile = "%s/trj%d.lh5"%(path,i)
    if os.path.exists(Trajfile):
        T = Trajectory.LoadFromLHDF(Trajfile)
        hct.CreateTrajFileWithHCstrings(T)
print "Done."
