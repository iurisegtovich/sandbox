import sys,os
import numpy as np
from scipy import savetxt,loadtxt
import msmbuilder
from msmbuilder import Trajectory,Serializer
 
def calculatedistance(AtomName1,ResidueID1,AtomName2,ResidueID2,trajfile,LongestTrajLength):
    """ Calculate the distance between given two atoms in given trajectory"""
    t = Trajectory.LoadFromLHDF(trajfile)
    Atom1 = (t['AtomNames'] == AtomName1)*(t['ResidueID'] == ResidueID1)
    Atom2 = (t['AtomNames'] == AtomName2)*(t['ResidueID'] == ResidueID2)
    distance = []
    for i in range(len(t['XYZList'])):
        x = (t['XYZList'][i,Atom1,:] - t['XYZList'][i,Atom2,:])[0]
        x = x.tolist()
        distance.append(np.dot(x,x)**0.5)
    distance += [-1]*(LongestTrajLength-len(t['XYZList']))
    return distance
	


#------------MAIN---------------


AtomName1 = 'C'
ResidueID1 = 1
AtomName2 = 'N'
ResidueID2 = 23
path ='/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster4.2'
Distances = []
ProjectInfo = Serializer.LoadFromHDF('%s/ProjectInfo.h5'%path)
LongestTrajLength = max(ProjectInfo['TrajLengths'])
os.chdir(path)
if os.path.exists('EndtoEndDistances.dat') :
    print "EndtoEndDistances.dat exists!"
    sys.exit()
print 'Calculating the eeDistance of each trajectory......' 
for i in range(ProjectInfo['NumTrajs']):
    trajfile = ProjectInfo['TrajFilePath']+ProjectInfo['TrajFileBaseName']+'%d'%i+ProjectInfo['TrajFileType']
    print '%d in %d Trajectories'%(i,ProjectInfo['NumTrajs']),trajfile
    d = calculatedistance(AtomName1,ResidueID1,AtomName2,ResidueID2,trajfile,LongestTrajLength)
    Distances.append(d)
print "Save data to ./EndtoEndDistance.dat"
savetxt('EndtoEndDistances.dat',Distances)
print "Done."

#distance = loadtxt('EndtoEndDistances.dat')
#print 'distance',distance
#distancemasked = np.ma.array(distance,mask=(distance==-1))
#print 'masked distance', distancemasked
#savetxt('EndtoEndDistances.Masked.dat',distancemasked)

