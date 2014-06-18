
import os,sys
import copy
import numpy as np
from scipy import savetxt,loadtxt
from msmbuilder import Trajectory


def FixGenFile(Mapping,GenFile,Outfile = './Gens.Fixed.lh5'):
    """
    Use Mapping.dat file get a fixed(mapped) generator file.
    New generator file will be Gens.Fixed.lh5
    """
    gen = Trajectory.LoadFromLHDF(GenFile)
    newgen = copy.deepcopy(gen)
    mapping = loadtxt(Mapping)
    GeneratorStateID = np.arange(len(gen['XYZList']))
    newgen['StateID'] = GeneratorStateID[mapping>=0]
    newgen['XYZList'] = gen['XYZList'][mapping>=0,:,:]   
    print "Save to %s"%Outfile
    newgen.SaveToLHDF(Outfile)

def GetSortedStateIDFromPopulations(PopulationFile='./Populations.dat',
                     NumOfPopStates=100,Order ='descend'):
    
    Populations = loadtxt(PopulationFile)
    PopulationsWithIndex = zip(Populations,range(len(Populations)))
    PopulationsWithIndex.sort()
    if Order == 'descend':
        PopulationsWithIndex.reverse()
    StateID = [PopulationsWithIndex[i][1] for i in range(NumOfPopStates)]
    
    return StateID

def ConvertMostPopulatedMappedStateStructuresToPDBFile(FixedGenFile,
                                                       PopulationFile='./Populations.dat',
                                                       NumOfPopStates=100,
                                                       Savepath = './generators'):
    
    StateID = GetSortedStateIDFromPopulations(PopulationFile,NumOfPopStates)
    genfixed = LoadTrajectory(FixedGenFile)
    newgen = copy.deepcopy(genfixed)
    newgen['XYZList'] = genfixed['XYZList'][StateID]
    CreatePath(Savepath)
    OutfileList = ['%s/gen%d.pdb'%(Savepath,i) for i in range(len(newgen['XYZList']))]
    pairfile = '%s/GenId-StateIdpair.txt'%Savepath
    print "Save the Generator Id and State Id to %s"%pairfile
    fn = open(pairfile,'w')
    fn.writelines("%s  State%d\n"%(OutfileList[i].split('/')[-1],StateID[i]) for i in range(len(OutfileList)))
    fn.close()
    ConvertToManyPDBFiles(newgen,OutfileList)
        
def CreatePath(Path):
    a = Path.split('/') 
    pathlist = ['/'.join(a[:i+1]) for i in range(len(a))]
    for path in pathlist:
        if not os.path.exists(path):
            print "Create directory:%s"%path
            os.mkdir(path)
            
    
def ConvertToManyPDBFiles(Trajectory,OutfileList):
    
    t = LoadTrajectory(Trajectory)
    if len(t['XYZList']) != len(OutfileList):
        raise("The output list do not match up the frames.")
    t1 = copy.deepcopy(t)
    t1['XYZList'] = -1*np.ones(t['XYZList'][:1].shape)
    for i in range(len(t['XYZList'])):
        t1['XYZList'][0,:,:] = t['XYZList'][i,:,:]
        print"Save to %s"%OutfileList[i]
        t1.SaveToPDB(OutfileList[i])
    print "Done."
    
def LoadTrajectory(trajectory):
    
    if isinstance(trajectory,str):
        try:
            t = Trajectory.LoadFromLHDF(trajectory)
            return t
        except IOError:
            raise IOError("Can not find %s"%trajectory)
    elif isinstance(trajectory,Trajectory):
        return trajectory

if __name__ == '__main__':
    Gens = '/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster4.2/Data/Gens.lh5'
    Mapping = '/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster4.2/lagtime50/Mapping.dat'
    #FixGenFile(Mapping,Gens)
    FixedGenFile = '/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster4.0/Data/Gens.lh5'
    Population = '/Users/tud51931/projects/surprisal/meansurprisalvalues-populationweighted-ERRRMSD4.0-700trajs.dat'
    ConvertMostPopulatedMappedStateStructuresToPDBFile(FixedGenFile,
                                                       PopulationFile=Population,
                                                       NumOfPopStates=20,
                                                       Savepath = './generators')        
        
