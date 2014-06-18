#calculate_observable_meanvalue
#export .dat file as an population input file of UseCKTest 

import os,sys
sys.path.append('/Users/tud51931/scripts/gfzhou')
from _predictobservables import PredictObservables
from scipy import savetxt


information = []
projectfile='/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster3.0/ProjectInfo.h5'
populationfile='/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster3.0/lagtime50/Populations.dat'
assignmentfile_fixed='/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster3.0/lagtime50/Assignments.Fixed.h5'
tmatrixfile='/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster3.0/lagtime50/tCounts.mtx'
rawdatafile='/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/result/RMSD.h5'


obj = PredictObservables(information,projectfile,populationfile,assignmentfile_fixed,tmatrixfile,rawdatafile)
obj.calculate_observable_meanvalue_state()
savetxt('AverageRMSD_States_RMSDcutoff3.0.dat',obj.ObservableMeanValue_State)
print "Save to AverageRMSD_States_RMSDcutoff3.0.dat"
