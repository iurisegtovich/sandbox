import os,sys
import numpy as np
from msmbuilder import Serializer
from scipy import savetxt,loadtxt
import matplotlib.pyplot as plt

#rmsd = loadtxt('/Users/tud51931/scripts/gfzhou/AverageRMSD_States_RMSDcutoff3.0.dat')
p= loadtxt('/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster3.0/lagtime50/Populations.dat')
sur = loadtxt('/Users/tud51931/scripts/gfzhou/SurprisalValue_State_100tcountsRMSD3.0.weighted.dat')
plt.figure()
plt.xlabel('StatePopulations')
plt.ylabel('WeightedSurprisalValue')
plt.xscale('log')
plt.yscale('log')
plt.ylim([0.0001,1])
plt.scatter(p,sur)
plt.title('WeightedSurprisalValue-StatePopulations')
plt.savefig('WeightedSurprisalValue-StatePopulations-RMSDcutoff3.0.100trajs.logscaleylim.png')
plt.show()

"""
#plt.hexbin(rmsd,sur,bins='log')
#plt.xlabel('RMSD(nm)')
#plt.ylabel('WeightedSurprisalValue')
#plt.title('WeightedSurprisalValue-RMSD')
#plt.savefig('WeightedSurprisalValue-RMSD-RMSDcutoff3.0.png')

plt.hexbin(p,sur,bins='log')
plt.xlabel('StatePopulations')
plt.ylabel('WeightedSurprisalValue')
plt.xlim([0,0.01])
plt.title('WeightedSurprisalValue-StatePopulations')
plt.savefig('WeightedSurprisalValue-StatePopulations-RMSDcutoff3.0.100trajs.xlim0.001.png')
plt.show()
"""


