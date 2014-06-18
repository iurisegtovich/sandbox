import os,sys
from scipy.io import mmread
from scipy import loadtxt
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
sys.path.append('/Users/tud51931/scripts/gfzhou')
import cutoffs as ct

def ChapmanKolmogorovTest(Tc_tau,Tc_ktau,Populations_tau,Populations_ktau,Mapping_tau,Mapping_ktau,k,cutoff,tau,states = 'all',mode=0):
    
    Probability_tau = 0
    Probability_ktau = 0
    if states == 'all' :
        P_tau = (np.matrix(Populations_tau)*(Tc_tau**k))[0]
        P_ktau = (np.matrix(Populations_ktau)*Tc_ktau)[0]
	
    else:
        Populations_tau_i = np.zeros(len(Populations_tau))
	Populations_ktau_i = np.zeros(len(Populations_ktau))
        for i in states:
	    if Mapping_tau[i] != -1 and Mapping_ktau[i] != -1:
		Populations_tau_i[i] = Populations_tau[int(Mapping_tau[i])]
		Populations_ktau_i[i] = Populations_ktau[int(Mapping_ktau[i])]
	    else:
		print "Warning:states %d does not match, remove state %d in the list" %(i,i)
		states.remove(i)
	Populations_tau_i /= Populations_tau_i.sum(0) #normalize
	Populations_ktau_i /=Populations_ktau_i.sum(0)
	P_tau = (np.matrix(Populations_tau_i)*(Tc_tau**k))[0]
        P_ktau = (np.matrix(Populations_ktau_i)*Tc_ktau)[0]
	print "Populations_ktau_i",Populations_ktau_i
	#print "P_tau,P_ktau", P_tau.shape,P_ktau
        for i in states:
	    Probability_tau += P_tau[0,int(Mapping_tau[i])]
	    Probability_ktau += P_ktau[0,int(Mapping_ktau[i])]
	print "return the Probability to stay in states"
	return Probability_tau,Probability_ktau 
				
            
               
        
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
    
#cutoff = ct.Cutoff(3,7,0.5)
#cutoffs = cutoff.cutoffs
cutoffs = [4.2]
Path = "/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/"
tau = 50
#stateslist = [0]
stateslist = list(np.linspace(0,200,201))
print "list",stateslist

for cutoff in cutoffs:
    RMSDdir = 'RMSDCluster%0.1f'%cutoff
    Tc_tau = mmread(os.path.join(Path,RMSDdir,'lagtime%d'%tau,'tProb.mtx'))
    Populations_tau = loadtxt(os.path.join(Path,RMSDdir,'lagtime%d'%tau,'Populations.dat'))
    Mapping_tau = loadtxt(os.path.join(Path,RMSDdir,'lagtime%d'%tau,'Mapping.dat'))
    Prob_tau,Prob_ktau = [1],[1]
    for k in range(1,9):
        ktau = k*tau
        Tc_ktau = mmread(os.path.join(Path,RMSDdir,'lagtime%d'%ktau,'tProb.mtx'))
        Populations_ktau = loadtxt(os.path.join(Path,RMSDdir,'lagtime%d'%ktau,'Populations.dat'))
        Mapping_ktau = loadtxt(os.path.join(Path,RMSDdir,'lagtime%d'%ktau,'Mapping.dat'))
        #E_ktau , E_ktau_rel = ChapmanKolmogorovTest(Tc_tau,Tc_ktau,Populations_tau,Populations_ktau,Mapping_tau,Mapping_ktau,k,cutoff,tau)
        Probability_tau, Probability_ktau = ChapmanKolmogorovTest(Tc_tau,Tc_ktau,Populations_tau,Populations_ktau,Mapping_tau,Mapping_ktau,k,cutoff,tau,states = stateslist)
        Prob_tau.append(Probability_tau)
        Prob_ktau.append(Probability_ktau)
    print "Probability_tau:",Prob_tau
    print "Probability_ktau;",Prob_ktau
    if 1:
	path = "/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/result/ChapmanKolmogorovTest"
        plt.figure()
        plt.plot(range(len(Prob_tau)),Prob_tau)
        plt.hold(True)
        plt.plot(range(len(Prob_ktau)),Prob_ktau)
        plt.legend(['Prob_tau','Prob_ktau'])
        plt.title("Cutoff = %0.1f tau = %d "%(cutoff,tau))
        
        if 0:
            plt.savefig("%s/tau%d_Cutoff%0.1f_k%d.png"%(path,tau,cutoff,k))        
        else:
            plt.show()
