import os,sys
import numpy as np
from msmbuilder import Serializer
from scipy import savetxt,loadtxt
import matplotlib.pyplot as plt


try:    
    R = Serializer.LoadFromHDF('/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster4.2/Data/RMSD.h5')
    rmsd = []
    for i in range(len(R['Data'])):
        for j in range(len(R['Data'][i])):
            if R['Data'][i,j] != -1:
                rmsd.append(R['Data'][i,j])
except IOError:
    print "Can't find RMSD.h5, please run CalculateProjectRMSD.py first to get RMSD.h5."
    sys.exit()
try:
    Nv = loadtxt('/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/result/numhelix_alltraj.txt')
    nv = []
    for i in range(len(Nv)):
        for j in range(len(Nv[i])):
            if Nv[i,j] != -1:
                nv.append(Nv[i,j])
except IOError:
    print "Can't find numhelix_alltraj.txt, please run computeNumhelix_alltrajs.py first."
    sys.exit()
    
plt.figure()
plt.hexbin(rmsd,nv,bins='log')
plt.xlabel('RMSD(nm)')
plt.ylabel('Number of Helix')
plt.title('Number of Helix - RMSD')
plt.savefig('Nv-RMSD.png')
#plt.show()

