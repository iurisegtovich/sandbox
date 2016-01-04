import os,sys
from scipy.io import mmread
from scipy import loadtxt
from msmbuilder import MSMLib
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

def ChapmanKolmogorovTest(Tc_tau,Tc_ktau,Populations_tau,Populations_ktau,Mapping_tau,Mapping_ktau,k,cutoff,tau,states = 'all'):
    
    if states=='all' :
        P_tau = (np.matrix(Populations_tau)*(Tc_tau**k))[0]
        P_ktau = (np.matrix(Populations_ktau)*Tc_ktau)[0]   
    else:
        np.zeros()
        for i in states:
    
    NumofStates = len(Mapping_tau)    
    E_ktau = np.zeros(NumofStates)
    E_ktau_rel = np.zeros(NumofStates)
   
    #print 'shape of P_tau', P_tau[0]
    for i in range(NumofStates):
        if Mapping_tau[i] != -1 and Mapping_ktau[i] != -1:     
           P1 = P_tau[0,int(Mapping_tau[i])]
           P2 = P_ktau[0,int(Mapping_ktau[i])]
           #print 'P1',P1
           #print 'P2',P2          
           E_ktau[i] = abs(P1-P2)
           E_ktau_rel[i] = abs(P1-P2)/P2
        else:
           E_ktau[i] = -1
           E_ktau_rel[i] = -1
           
    if 1:
        path = "/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/result/ChapmanKolmogorovTest"
        plt.figure()
        plt.plot(range(NumofStates),E_ktau)
        plt.hold(True)
        plt.plot(range(NumofStates),E_ktau_rel)
        plt.legend(['E_%dtau'%k,'E_%dtau_rel'%k])
        plt.title("Cutoff = %0.1f tau = %d k = %d"%(cutoff,tau,k))
        plt.xaxil
        if 1:
            plt.savefig("%s/tau%d_Cutoff%0.1f_k%d.png"%(path,tau,cutoff,k))        
        else:
            plt.show()
        
    return E_ktau , E_ktau_rel
    
__metaclass__= type 

class Cutoff:
    
    def init(self,startcutoff,endcutoff,stepsize):
        self.cutoffs = []
        i = startcutoff
        while i <= endcutoff:
            self.cutoffs.append(round(i,1))
            i += stepsize
        return self.cutoffs
    
    def ext(self,startcutoff,endcutoff,stepsize):
        newcutoffs = np.arange(startcutoff,endcutoff,stepsize)
        oldcutoffs = self.cutoffs[:]
        for x in newcutoffs:
            if not self.cutoffs.count(x):
                self.cutoffs.append(round(x,1))
        self.cutoffs.sort()
        return self.cutoffs



cutoff = Cutoff()
cutoffs = cutoff.init(4.1,4.5,0.1) #(exclude stop point that is 8.5)
#cutoffs = cutoff.ext(2,2.5,0.1)
print "cutoffs",cutoffs

Path = "/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/"
tau = 100

for cutoff in cutoffs:
    RMSDdir = 'RMSDCluster%0.1f'%cutoff
    Tc_tau = mmread(os.path.join(Path,RMSDdir,'lagtime%d'%tau,'tProb.mtx'))
    Populations_tau = loadtxt(os.path.join(Path,RMSDdir,'lagtime%d'%tau,'Populations.dat'))
    Mapping_tau = loadtxt(os.path.join(Path,RMSDdir,'lagtime%d'%tau,'Mapping.dat'))
    for k in range(2,6):
        ktau = k*tau
        Tc_ktau = mmread(os.path.join(Path,RMSDdir,'lagtime%d'%ktau,'tProb.mtx'))
        Populations_ktau = loadtxt(os.path.join(Path,RMSDdir,'lagtime%d'%ktau,'Populations.dat'))
        Mapping_ktau = loadtxt(os.path.join(Path,RMSDdir,'lagtime%d'%ktau,'Mapping.dat'))
        E_ktau , E_ktau_rel = ChapmanKolmogorovTest(Tc_tau,Tc_ktau,Populations_tau,Populations_ktau,Mapping_tau,Mapping_ktau,k,cutoff,tau)
        print "k,E_ktau_rel",k,E_ktau_rel
    