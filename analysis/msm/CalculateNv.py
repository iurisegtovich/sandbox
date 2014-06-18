
import os,sys
import argparse
import numpy as np
from scipy import savetxt
from msmbuilder.geometry import dihedral
from msmbuilder import Trajectory,Serializer

sys.path.append('/Users/tud51931/scripts/gfzhou')
import HelixCoilTools as hct



parser = argparse.ArgumentParser()
parser.add_argument('project',help="Path to ProjectInfo.h5,default=ProjectInfo.h5",default="ProjectInfo.h5")
parser.add_argument('-o','--Output',help="Output file. default=Nv.dat",default="Nv.dat")
args = parser.parse_args()


if os.path.exists('Nv.dat') :
    print "Nv.dat exists!"
    sys.exit()


ProjectInfo = Serializer.LoadFromHDF(args.project)
LongestTrajLength = max(ProjectInfo['TrajLengths'])
NumberOfHelix = -1*np.ones((ProjectInfo['NumTrajs'],LongestTrajLength))
print 'Calculating the Number of Helix for each trajectory......' 
for i in range(ProjectInfo['NumTrajs']): 
    trajfile = ProjectInfo['TrajFilePath']+ProjectInfo['TrajFileBaseName']+'%d'%i+ProjectInfo['TrajFileType']
    if os.path.exists(trajfile):
        print '%d in %d Trajectories'%(i,ProjectInfo['NumTrajs']),trajfile
        t = Trajectory.LoadFromLHDF(trajfile)
        Nv = hct.compute_numhelix_trajectory(t)
        NumberOfHelix[i,:len(Nv)] = Nv[:]
print "Save to %s"%args.Output
savetxt(args.Output,NumberOfHelix)
print "Done."