import os,sys
from random import choice
import numpy as np
from msmbuilder import Serializer,Trajectory
from scipy import savetxt,loadtxt
from scipy.io import mmread
import matplotlib.pyplot as plt
import pickle

sys.path.append("/Users/tud51931/scripts/gfzhou/")
import HelixCoilTools as hct
import CalculateRg


def Nvprediction():
    try:
        mean_numhelix_states = loadtxt('mean_numhelix_states.dat')
    except IOError:
    
        StatesAsi = hct.get_StatesAssignments(Assignments)
        NumHelix_states = hct.compute_numhelix_states(StatesAsi)
        #print "NumHelix_states",NumHelix_states
        #savetxt('NumHelix_states',NumHelix_states)
        states = [int(i) for i in NumHelix_states.keys()]
        states.sort()
        mean_numhelix_states = []
        std_numhelix_states = []
        for state in states:
            mean_numhelix_states.append(np.mean(NumHelix_states['%d'%state]))
            std_numhelix_states.append(np.std(NumHelix_states['%d'%state]))
        savetxt('mean_numhelix_states.dat',mean_numhelix_states)
        savetxt('std_numhelix_states.dat',std_numhelix_states)
    #plt.figure()
    #plt.errorbar(states,mean_numhelix_states,std_numhelix_states)
    #plt.xlabel("State ID")
    #plt.ylabel("Number of Helix")
    #plt.savefig("Numhelix_states")
    #plt.show()
    
    P0 = np.zeros(len(Population))
    for data in Assignments['Data']:
        P0[data[0]] += 1
    P0 = P0/P0.sum()
    populationslist = []
    for k in range(140):     # tau = 50, so 140*50 = 7000
        populationslist.append(P0)
        P0 *= Tmatrix
     
    numhelix = np.dot(np.array(populationslist), np.array(mean_numhelix_states).reshape(-1,1))
    print numhelix
    numhelix = numhelix.reshape(1,-1)[0]
    plt.figure()
    plt.plot(np.arange(0,7000,50),numhelix,'ro',)# tau = 50, so 140*50 = 7000
    plt.hold(True)
    
    Counts = -1*np.ones((ProjectInfo['NumTrajs'],max(ProjectInfo['TrajLengths'])))
    print Counts.shape
    
    for i in range(0,93):
        T = Trajectory.LoadFromHDF('/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/sourcedata/Trajectories/trj%d_hc.h5'%i)
        Hcount = hct.count_Helix(T)
        Counts[i,:len(Hcount)] = Hcount[:]
    
    Counts_ma = np.ma.array(Counts,mask=[Counts==-1])
    H_mean = Counts_ma.mean(0)
    H_std = Counts_ma.std(0)
    print H_mean
    
    plt.plot(range(len(H_mean)),H_mean,'b')
    
    plt.title('Nv-steps')
    plt.xlabel('Steps')
    plt.ylabel('Nv')
    plt.legend(('Nv_msm','Nv_rawdata'),loc = 'upper left')
    figname = 'Nv_prediction_%sCluster%0.1f_tau%d.png'%(metrics,cutoff,tau)
    plt.savefig(figname)
    print "Save to %s"%figname
    #plt.show()
    
def Rgprediction():
    try:
        Rgs = loadtxt('/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/result/Rgs.dat')
    except IOError:
        print "Can't find Rgs.dat, please run CalculateRg.py first."
        sys.exit()
    StatesAsi = hct.get_StatesAssignments(Assignments)
    Rgs_states = {}
    for state in StatesAsi.keys():
        for trajid in StatesAsi[state].keys():
            for frame in StatesAsi[state][trajid]:
                Rgs_states.setdefault(state,[]).append(Rgs[int(trajid)][int(frame)])
            
    states = [int(i) for i in Rgs_states.keys()]
    states.sort()
    mean_rg_states = []
    std_rg_states = []
    for state in states:
        mean_rg_states.append(np.mean(Rgs_states['%d'%state]))
        std_rg_states.append(np.std(Rgs_states['%d'%state]))
    #savetxt('mean_numhelix_states0.dat',mean_numhelix_states)
    #savetxt('std_numhelix_states0.dat',std_numhelix_states)
    print mean_rg_states
    
    P0 = np.zeros(len(Population))
    for data in Assignments['Data']:
        P0[data[0]] += 1
    P0 = P0/P0.sum()
    populationslist = []
    for k in range(140):
        populationslist.append(P0)
        P0 *= Tmatrix    
        
    Rgs_predicted = np.dot(np.array(populationslist), np.array(mean_rg_states).reshape(-1,1))
    print Rgs_predicted
    Rgs_predicted = Rgs_predicted.reshape(1,-1)[0]
    plt.figure()
    plt.plot(np.arange(0,7000,50),Rgs_predicted,'ro',)
    plt.hold(True)
    
    
    Counts_ma = np.ma.array(Rgs,mask=[Rgs==-1])
    Rgs_mean = Counts_ma.mean(0)
    Rgs_std = Counts_ma.std(0)
    print Rgs_mean
    
    plt.plot(range(len(Rgs_mean)),Rgs_mean,'b')
    
    plt.title('Rgs-steps')
    plt.xlabel('Steps')
    plt.ylabel('Rgs')
    plt.legend(('Rgs_msm','Rgs_rawdata'),loc = 'upper left')
    figname = 'Rgs_prediction_%sCluster%0.1f_tau%d.png'%(metrics,cutoff,tau)
    plt.savefig(figname)
    print "Save to N%s"%figname
    #plt.show()    

def RMSDprediction():
    try:    
        R = Serializer.LoadFromHDF('/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster4.2/Data/RMSD.h5')
    except IOError:
        print "Can't find RMSD.h5, please run CalculateProjectRMSD.py first to get RMSD.h5."
        sys.exit()    
    RMSD = R['Data']
    StatesAsi = hct.get_StatesAssignments(Assignments)
    RMSD_states = {}
    for state in StatesAsi.keys():
        for trajid in StatesAsi[state].keys():
            for frame in StatesAsi[state][trajid]:
                RMSD_states.setdefault(state,[]).append(RMSD[int(trajid)][int(frame)])    
                
    states = [int(i) for i in RMSD_states.keys()]
    states.sort()
    mean_rmsd_states = []
    std_rmsd_states = []
    for state in states:
        mean_rmsd_states.append(np.mean(RMSD_states['%d'%state]))
        std_rmsd_states.append(np.std(RMSD_states['%d'%state]))
    #savetxt('mean_numhelix_states0.dat',mean_numhelix_states)
    #savetxt('std_numhelix_states0.dat',std_numhelix_states)
    print mean_rmsd_states
    
    P0 = np.zeros(len(Population))
    for data in Assignments['Data']:
        P0[data[0]] += 1
    P0 = P0/P0.sum()
    populationslist = []
    for k in range(140):
        populationslist.append(P0)
        P0 *= Tmatrix    
        
    RMSD_predicted = np.dot(np.array(populationslist), np.array(mean_rmsd_states).reshape(-1,1))
    print RMSD_predicted
    RMSD_predicted = RMSD_predicted.reshape(1,-1)[0]
    plt.figure()
    plt.plot(np.arange(0,7000,50),RMSD_predicted,'ro',)
    plt.hold(True)
    
    
    Counts_ma = np.ma.array(RMSD,mask=[RMSD==-1])
    RMSD_mean = Counts_ma.mean(0)
    RMSD_std = Counts_ma.std(0)
    print RMSD_mean
    
    plt.plot(range(len(RMSD_mean)),RMSD_mean,'b')
    
    plt.title('RMSD-steps')
    plt.xlabel('Steps')
    plt.ylabel('RMSD')
    plt.legend(('RMSD_msm','RMSD_rawdata'),loc = 'upper right')
    figname = 'RMSD_prediction_%sCluster%0.1f_tau%d.png'%(metrics,cutoff,tau)
    plt.savefig(figname)
    print "Save to %s"%figname
    #plt.show()        
    
def EEdistanceprediction():
    try:    
        EEd = loadtxt('/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/result/EndtoEndDistances.dat')
    except IOError:
        print "Can't find EndtoEndDistances.dat, please run CalculateEndtoEndDistance.py first to get EndtoEndDistances.dat."
        sys.exit()    
    StatesAsi = hct.get_StatesAssignments(Assignments)
    EEd_states = {}
    for state in StatesAsi.keys():
        for trajid in StatesAsi[state].keys():
            for frame in StatesAsi[state][trajid]:
                EEd_states.setdefault(state,[]).append(EEd[int(trajid)][int(frame)])    
                
    states = [int(i) for i in EEd_states.keys()]
    states.sort()
    mean_eed_states = []
    std_eed_states = []
    for state in states:
        mean_eed_states.append(np.mean(EEd_states['%d'%state]))
        std_eed_states.append(np.std(EEd_states['%d'%state]))
    #savetxt('mean_numhelix_states0.dat',mean_numhelix_states)
    #savetxt('std_numhelix_states0.dat',std_numhelix_states)
    print mean_eed_states
    
    P0 = np.zeros(len(Population))
    for data in Assignments['Data']:
        P0[data[0]] += 1
    P0 = P0/P0.sum()
    populationslist = []
    for k in range(140):
        populationslist.append(P0)
        P0 *= Tmatrix    
        
    EEd_predicted = np.dot(np.array(populationslist), np.array(mean_eed_states).reshape(-1,1))
    print EEd_predicted
    EEd_predicted = EEd_predicted.reshape(1,-1)[0]
    plt.figure()
    plt.plot(np.arange(0,7000,50),EEd_predicted,'ro',)
    plt.hold(True)
    
    
    Counts_ma = np.ma.array(EEd,mask=[EEd==-1])
    EEd_mean = Counts_ma.mean(0)
    EEd_std = Counts_ma.std(0)
    print EEd_mean
    
    plt.plot(range(len(EEd_mean)),EEd_mean,'b')
    
    plt.title('EEd-steps')
    plt.xlabel('Steps')
    plt.ylabel('EEd')
    plt.legend(('EEd_msm','EEd_rawdata'),loc = 'upper left')
    figname = 'EEd_prediction_%sCluster%0.1f_tau%d.png'%(metrics,cutoff,tau)
    plt.savefig(figname)
    print "Save to %s"%figname
    #plt.show()            
    
def create_hcstrings_states(Assignments,outfile = 'HCstrings_states.txt'):
    SA = hct.get_StatesAssignments(Assignments)
    states = SA.keys()
    HCstrings_states = {}
    n = 0
    for state in states:
        n +=1
        print "Get HC strings for state %d/%d"%(n,len(states))
        TrajID = SA[state].keys()
        numhelix_state = []
        HCstrings_states[state] = []
        for trajid in TrajID:          
            TrajFile = '/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/sourcedata/Trajectories/trj%s_hc.lh5'%trajid
            Traj = Trajectory.LoadFromLHDF(TrajFile)
            HCstrings_states[state] += [Traj['HCs'][i] for i in SA[state][trajid]] 
    fn = outfile    
    if os.path.exists(fn):
        newfn = fn + '.bck'
        os.system('mv %s %s'%(fn,newfn))
    print "Write HCstings of states into %s"%fn        
    HCfile = open(fn,'w')
    pickle.dump(HCstrings_states,HCfile)
    HCfile.close()
    print "Done."
   
def calculate_SequenceEntropy(SequenceList,ID = 'h'):
    S = SequenceList
    NumSequences = len(S)
    SequenceLength = len(S[0])
    Populations_positions = []
    for j in range(SequenceLength):
        NumIDs = 0.0        
        for i in range(NumSequences):
            if S[i][j] == ID:
                NumIDs += 1.0
        Populations_positions.append(NumIDs/NumSequences)
    SequenceEntropy = np.sum([-Populations_positions[k]*np.log(Populations_positions[k]) for k in range(SequenceLength)])
    return SequenceEntropy

def SequenceEntropy_states(inputfile='HCstrings_states.txt'):
    fn = open(inputfile,'r')
    HCstrings_states = pickle.load(fn)
    fn.close()
    states = [int(i) for i in HCstrings_states.keys()]
    StatesEntropy = -1*np.ones(max(states)+1)
    for i in states:
        print "Calculating %d of %d states"%(i,len(states)) 
        StatesEntropy[i] = calculate_SequenceEntropy(HCstrings_states[str(i)])
    
    return StatesEntropy

def barchartsforStatesEntropy():
    StatesEntropy = SequenceEntropy_states()
    N = len(StatesEntropy)
    ind = np.arange(N)*2  # the x locations for the groups
    width = 0.35          
    plt.figure()
    rects1 = plt.bar(ind, StatesEntropy, width, color='r')
    plt.ylabel('Sequence Entropy')
    plt.title('Sequence Entropy over States')
    
    plt.xticks(ind+width, ['S%d'%i for i in range(N)] )
    
    #autolabel(rects1)
       
    plt.show()        
    
def autolabel(rects):
    # attach some text labels
    for rect in rects:
        height = rect.get_height()
        plt.text(rect.get_x()+rect.get_width()/2., 1.05*height, '%d'%int(height),
                ha='center', va='bottom')

def GetRgsforGeneratorFile():
    Output = '%s/Rgs-gens.dat'%os.getcwd()
    Trajfiles = []
    trjfile = Gens
    Trajfiles.append(trjfile)
    CalculateRg.CalculateRg(Trajfiles,Output)

def GetHCStringsforTrajectory(trajectory):
    if isinstance(trajectory,str):
        gens = Trajectory.LoadFromLHDF(trajectory)
    dihedrals = hct.ComputeDihedralsFromTrajectory(gens)
    HCs = hct.ConvertDihedralsToHCStrings(dihedrals)
    print HCs

    
#def GetHCStringsforProject(ProjectInfo)

if __name__ == '__main__':
    tau = 50
    cutoff = 4.2
    metrics = 'rmsd'
    
    if metrics.lower() == 'dihedral':
        Path = "/Users/tud51931/projects/MSM/msm/ff03-dihedralhybrid"
        metrics = 'Dihedral'
    elif metrics.lower() == 'rmsd':
        Path = "/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter"
        metrics = 'RMSD'    
    Path = os.path.join(Path,'%sCluster%0.1f'%(metrics,cutoff))
        
    ProjectInfo = Serializer.LoadFromHDF('%s/ProjectInfo.h5'%Path)
    Population = loadtxt('%s/lagtime%d/Populations.dat'%(Path,tau))
    Assignments = Serializer.LoadFromHDF("%s/lagtime%d/Assignments.Fixed.h5"%(Path,tau))
    Tmatrix = mmread('%s/lagtime%d/tProb.mtx'%(Path,tau))
    Gens = '%s/Data/Gens.lh5'%Path    

    Nvprediction()
    EEdistanceprediction()
    RMSDprediction()
    Rgprediction()    
    
    
    #barchartsforStatesEntropy()
    #SequenceEntropy_states()
    #GetRgsforGeneratorFile()
    #GetHCStringsforTrajectory(Gens)
