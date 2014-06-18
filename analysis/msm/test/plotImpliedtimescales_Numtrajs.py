import os,sys
from scipy import loadtxt
import numpy as np
import matplotlib.pyplot as plt

Impliedtime_mean = []
Impliedtime_std = []
for i in range(2,87):
    path = '/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster4.2/lagtime50/bootstrap/%dtrajs'%i
    time = []    
    for j in range(1,101):
        fn = os.path.join(path,'ImpliedTimescales_%d.dat'%j)
        time.append(loadtxt(fn)[0])
    time_mean = np.mean(time)
    time_std = np.std(time)
    time_sam_std = time_std/(len(time)**0.5)
    Impliedtime_mean.append(time_mean)
    Impliedtime_std.append(time_sam_std)

plt.figure()
plt.errorbar(range(2,87),Impliedtime_mean,Impliedtime_std)
plt.xlabel('Number of Trajectories')
plt.ylabel('Impliedtimescales')
plt.savefig('/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster4.2/result/Impliedtimes_Numtrajs.png')
plt.show()
