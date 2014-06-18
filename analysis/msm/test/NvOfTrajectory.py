import os,sys
import numpy as np
from msmbuilder import Serializer,Trajectory
import matplotlib.pyplot as plt
sys.path.append("/Users/tud51931/scripts/gfzhou")
import HelixCoilTools as hct

ProjectInfo = Serializer.LoadFromHDF('/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/sourcedata/ProjectInfo.h5')
Counts = -1*np.ones((ProjectInfo['NumTrajs'],max(ProjectInfo['TrajLengths'])))
print Counts.shape

Savepath = '/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/result/NvOfTrajectory'
plt.figure()
plt.xlabel('Steps')
plt.ylabel('Nv')
plt.hold(False)
for i in range(0,93):
    T = Trajectory.LoadFromHDF('/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/sourcedata/Trajectories/trj%d_hc.h5'%i)
    Hcount = hct.count_Helix(T)
    plt.title('Nv-steps of Traj%d'%i)
    plt.plot(range(len(Hcount)),Hcount,'.')
    print 'Save figure to %s/Nvoftraj%d.png'%(Savepath,i)
    plt.savefig('%s/Nvoftraj%d.png'%(Savepath,i))
    Counts[i,:len(Hcount)] = Hcount[:]

Counts_ma = np.ma.array(Counts,mask=[Counts==-1])
H_mean = Counts_ma.mean(0)
H_std = Counts_ma.std(0)
print H_mean

plt.figure()
plt.plot(range(len(H_mean)),H_mean,'b')
plt.title('AverageNv-Steps Of All Trajectories')
plt.xlabel('Steps')
plt.ylabel('AverageNv')
print 'Save figure to %s/NvofAlltrajs.png'%Savepath
plt.savefig('%s/NvofAlltrajs.png'%Savepath)

