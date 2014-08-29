import os,sys
import argparse
from mdtraj import io
from scipy import loadtxt,savetxt
import matplotlib.pyplot as plt

def GetSortedStateIDFromPopulations(PopulationFile='./Populations.dat',
                     NumOfPopStates=100,Order ='descend'):
    
    Populations = loadtxt(PopulationFile)
    PopulationsWithIndex = zip(Populations,range(len(Populations)))
    PopulationsWithIndex.sort()
    if Order == 'descend':
        PopulationsWithIndex.reverse()
    StateID = [PopulationsWithIndex[i][1] for i in range(NumOfPopStates)]
    
    return StateID
     
def locateGensOnPlot(PopulationFile='./Populations.dat',
                     NumOfPopStates=100,Order ='descend',
                     StateID=None):
    
    try:
        StateID_Populations = GetSortedStateIDFromPopulations(PopulationFile,NumOfPopStates,Order)
        PlotBasedOnStateID(StateID_Populations)
    except IOError:
        if StateID == None:
            print "Can't find Population file, and StateID is None.Please offer at least one of them."%PopulationFile
            sys.exit()
    PlotBasedOnStateID(StateID_Populations)      
    
def PlotBasedOnStateID(StateID=None):     
    
    try:    
        R = loadtxt('/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster4.0/RMSD-gens.dat')
        rmsd_gens = []
        for i in range(len(R)):
            if R[i] != -1:
                rmsd_gens.append(R[i])
        rmsd_gens = np.array(rmsd_gens)
    except IOError:
        print "Can't find RMSD-gens.dat, please run CalculateRMSD.py first to get RMSD-gens.dat."
        raise IOError
    try:
        Rgs = loadtxt('/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster4.0/Rgs-gens.dat')
        rgs_gens = []
        for i in range(len(Rgs)):
            if Rgs[i] != -1:
                rgs_gens.append(Rgs[i])
        rgs_gens = np.array(rgs_gens)    
    except IOError:
        print "Can't find Rgs-gens.dat, please run CalculateRg.py first."
        raise IOError    

    if StateID is not None:
        plt.hold(True) 
        for i in StateID:
            plt.plot(rmsd_gens[i],rgs_gens[i],'ro')
            plt.text(rmsd_gens[i],rgs_gens[i],'%d'%i,fontsize=8)
            
parser = argparse.ArgumentParser()
parser.add_argument('-rmsd','--RMSD',help="Input RMSD.h5 file",required=True)
parser.add_argument('-rg','--Rg',help="Input Rgs.dat file",required=True)
parser.add_argument('-l','--Locate',help="Locate states on the Rg-RMSD graph",type=str)
parser.add_argument('-o','--Output',help="Output file (graph) name.Default: Rg-RMSD.png",default="Rg-RMSD.png")

args = parser.parse_args()

try:    
    R = io.loadh(args.RMSD,'arr_0')
    rmsd = []
    for i in range(R.shape[0]):
        for j in range(R.shape[1]):
            if R[i,j] != -1:
                rmsd.append(R[i,j])
except IOError:
    print "Can't find RMSD.h5"
    raise IOError
try:
    Rgs = io.loadh(args.Rg,'arr_0')
    rgs = []
    for i in range(Rgs.shape[0]):
        for j in range(Rgs.shape[1]):
            if Rgs[i,j] != -1:
                rgs.append(Rgs[i,j])
except IOError:
    print "Can't find Rgs.h5."
    raise IOError

plt.figure()
plt.hexbin(rmsd,rgs,bins='log')
plt.xlabel('RMSD(nm)')
plt.ylabel('Rg(nm)')
plt.title('Rg-RMSD')
if args.Locate !=None:
    locateGensOnPlot(PopulationFile=args.Locate,
                 NumOfPopStates=20,Order ='descend',
                 StateID=[0])
plt.savefig(args.Output)
print "Wrote: %s"%args.Output
#plt.show()
