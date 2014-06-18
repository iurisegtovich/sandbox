import sys,os
sys.path.append('Users/tud51931/scripts/gfzhou')
import ExtractStateStructures as ES
import ChapmanKolmogorovTestBootstrap as CKTest

#cutoffs = [4.0,4.1,4.2,4.3,4.4,4.5]
cutoffs = [3.0]
tau = 50
"""
for cutoff in cutoffs:
    RMSDdir='/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster%0.1f'%cutoff
    Populationfile = os.path.join(RMSDdir,'lagtime%d'%tau,'Populations.dat')
    StateID = ES.GetSortedStateIDFromPopulations(Populationfile,NumOfPopStates=100,Order ='descend')
    CKTest.ChapmanKolmogorovTest(RMSDdir,stateslist=StateID,tau=tau,kmax=8,cutoff=cutoff,bootstrap=False)
"""
for cutoff in cutoffs:
    RMSDdir='/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster%0.1f'%cutoff
    Populationfile = ('/Users/tud51931/scripts/gfzhou/AverageRMSD_States_RMSDcutoff3.0.dat')
    StateID = ES.GetSortedStateIDFromPopulations(Populationfile,NumOfPopStates=400,Order ='ascend')
    CKTest.ChapmanKolmogorovTest(RMSDdir,stateslist=StateID,tau=tau,kmax=8,cutoff=cutoff,bootstrap=False)