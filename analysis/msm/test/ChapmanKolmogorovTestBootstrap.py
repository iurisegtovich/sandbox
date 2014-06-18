"""
Chapman-Kolmogrov Test is a method to test and validate a Markov State Model.
Reference: 
1. Martin Senne, Benjamin Trendelkamp-Schroer, Antonia S.J.S. Mey, Christof Schutte, and Frank Noe,"EMMA: A Software Package for Markov Model Building and Analysis",J. Chem. Theory Comput. 2012, 8, 2223-2238.

----------
CHANGE LOG
08,OCT,2012 Guangfeng Zhou
----------

Script Author:
Guangfeng Zhou
Graduate Student
Dr.Voelz Lab
Chemistry Department
Temple University

"""

import os,sys
from scipy.io import mmread
from scipy import loadtxt
import numpy as np
import matplotlib.pyplot as plt


"""
Chapman-Kolmogorov Test
Input: 
    T_tau,T_ktau,Equilibrium Populations_tau,Equilibrium Populations_ktau_,Mapping_tau,Mapping_ktau,k,tau,selected states
    
    Selected states: A list contains the number of selected states.If some states are not used in both Transition Probability matrices(T_tau and T_ktau), they will be removed and a fixed list will be generate which you can print it out to see the states this test actually used.
    
    T--tProb.mtx
    Equilibrium Populations -- Populations.dat
    Mapping -- Mapping.dat
    
    

Output: 
    Probability_tau -- Probability to stay at selected states after time k*tau, predicted from MSMs P(tau)*(T(tau)**k)
    Probability_ktau -- Probability to stay at selected states after time k*tau calculated from raw data P(ktau)*T(ktau)

"""


def CalculateStatesProbability(Tc_tau,Tc_ktau,Populations_tau,Populations_ktau,Mapping_tau,Mapping_ktau,k,cutoff,tau,states):
    
    Probability_tau = 0
    Probability_ktau = 0
    Populations_tau_i = np.zeros(len(Populations_tau))
    Populations_ktau_i = np.zeros(len(Populations_ktau))
    for i in states:
	try: 
	    if Mapping_tau[i] != -1 and Mapping_ktau[i] != -1:
		Populations_tau_i[i] = Populations_tau[int(Mapping_tau[i])]
		Populations_ktau_i[i] = Populations_ktau[int(Mapping_ktau[i])]
	    else:
		print "Warning:states %d does not match, remove state %d in the list" %(i,i)
		states.remove(i)
	except IndexError:
	    print "state number %d is out of bound, remove state %d in the list"%(i,i)
	    states.remove(i)
    #Normalize the initial populations
    Populations_tau_i /= Populations_tau_i.sum(0) 
    Populations_ktau_i /=Populations_ktau_i.sum(0)
    #P_tau = (np.matrix(Populations_tau_i)*(Tc_tau**k))[0]
    #P_ktau = (np.matrix(Populations_ktau_i)*Tc_ktau)[0]
    P_tau = (np.matrix(Populations_tau_i)*(Tc_tau**k))[0]
    P_ktau = (np.matrix(Populations_ktau_i)*Tc_ktau)[0]
    for i in states:
	Probability_tau += P_tau[0,int(Mapping_tau[i])]
	Probability_ktau += P_ktau[0,int(Mapping_ktau[i])]
	
    return Probability_tau,Probability_ktau 
				
            
def ensemble_mean_error(Tc_tau,Tc_ktau,Populations_tau,Populations_ktau,Mapping_tau,Mapping_ktau,k,cutoff,tau):
    """
    TODO: rewrite the scripts, P_tau should be the prediction from model, P_ktau should be raw data which is just Populations_ktau.
    
    the function name should be RMSEinStatePopulations(TJ paper)
    """
    NumofStates = len(Mapping_tau)    
    E_ktau = np.zeros(NumofStates)
    E_ktau_rel = np.zeros(NumofStates)
    P_tau = (np.matrix(Populations_tau)*(Tc_tau**k))[0]
    P_ktau = (np.matrix(Populations_ktau)*Tc_ktau)[0]    
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
           
    if 0:
        path = "/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/result/ChapmanKolmogorovTest"
        plt.figure()
        plt.plot(range(NumofStates),E_ktau)
        plt.hold(True)
        plt.plot(range(NumofStates),E_ktau_rel)
        plt.legend(['E_%dtau'%k,'E_%dtau_rel'%k])
        plt.title("Cutoff = %0.1f tau = %d k = %d"%(cutoff,tau,k))
        plt.ylim(0,1)
        if 1:
            plt.savefig("%s/tau%d_Cutoff%0.1f_k%d.png"%(path,tau,cutoff,k))        
        else:
            plt.show()
        
    return E_ktau , E_ktau_rel
    
def ChapmanKolmogorovTest(RMSDdir='/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster4.2',
                          stateslist=[],tau=50,kmax=8,cutoff=4.2,bootstrap=False,bootstrapnumber=100):
    
    Prob_tau,Prob_ktau = [1],[1]
    Prob_tau_std,Prob_ktau_std = [0],[0]
    for k in range(1,kmax+1):
	if bootstrap:
	    ktau = k*tau
	    filepath_tau = os.path.join(RMSDdir,'lagtime%d'%tau,'bootstrap')
	    filepath_ktau = os.path.join(RMSDdir,'lagtime%d'%ktau,'bootstrap')    
	    Prob_tau_bootstrap = []
	    Prob_ktau_bootstrap = []
	
	    for i in range(1,bootstrapnumber+1):
		print "k = %d of 8 , bootstrap directory = %d of 100"%(k,i)
		Tc_tau = mmread(os.path.join(filepath_tau,'%d'%i,'tProb.mtx'))
		Populations_tau = loadtxt(os.path.join(filepath_tau,'%d'%i,'Populations.dat'))
		Mapping_tau = loadtxt(os.path.join(filepath_tau,'%d'%i,'Mapping.dat'))	
		
		Tc_ktau = mmread(os.path.join(filepath_ktau,'%d'%i,'tProb.mtx'))
		Populations_ktau = loadtxt(os.path.join(filepath_ktau,'%d'%i,'Populations.dat'))
		Mapping_ktau = loadtxt(os.path.join(filepath_ktau,'%d'%i,'Mapping.dat'))
		
		Probability_tau, Probability_ktau = CalculateStatesProbability(Tc_tau,Tc_ktau,Populations_tau,Populations_ktau,Mapping_tau,Mapping_ktau,k,cutoff,tau,states = stateslist)
		Prob_tau_bootstrap.append(Probability_tau)
		Prob_ktau_bootstrap.append(Probability_ktau)
	    Prob_tau.append(np.mean(Prob_tau_bootstrap))
	    Prob_tau_std.append(np.std(Prob_tau_bootstrap))
	    Prob_ktau.append(np.mean(Prob_ktau_bootstrap))
	    Prob_ktau_std.append(np.std(Prob_ktau_bootstrap))
	    
	    
	else:
	    ktau = k*tau
	    filepath_tau = os.path.join(RMSDdir,'lagtime%d'%tau)
	    filepath_ktau = os.path.join(RMSDdir,'lagtime%d'%ktau)    
	    Prob_tau_bootstrap = []
	    Prob_ktau_bootstrap = []
	
	    Tc_tau = mmread(os.path.join(filepath_tau,'tProb.mtx'))
	    Populations_tau = loadtxt(os.path.join(filepath_tau,'Populations.dat'))
	    Mapping_tau = loadtxt(os.path.join(filepath_tau,'Mapping.dat'))	
	    
	    Tc_ktau = mmread(os.path.join(filepath_ktau,'tProb.mtx'))
	    Populations_ktau = loadtxt(os.path.join(filepath_ktau,'Populations.dat'))
	    Mapping_ktau = loadtxt(os.path.join(filepath_ktau,'Mapping.dat'))
	    
	    Probability_tau, Probability_ktau = CalculateStatesProbability(Tc_tau,Tc_ktau,Populations_tau,Populations_ktau,Mapping_tau,Mapping_ktau,k,cutoff,tau,states = stateslist)
	    
	    
	    Prob_tau.append(Probability_tau)
	    Prob_ktau.append(Probability_ktau)    
	    
	    
	print "list-fixed",stateslist

    if bootstrap:
	if 1:
	    path = "/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/result/ChapmanKolmogorovTest"
	    plt.figure()
	    plt.errorbar(range(len(Prob_tau)),Prob_tau,Prob_tau_std)
	    plt.hold(True)
	    plt.errorbar(range(len(Prob_ktau)),Prob_ktau,Prob_ktau_std)
	    plt.legend(['Prob_tau','Prob_ktau'])
	    plt.title("Cutoff = %0.1f tau = %d "%(cutoff,tau))
	    plt.ylim(0,1)
	    if 1:
		plt.savefig("CKtestbootstrap_tau%d_Cutoff%0.1f.png"%(tau,cutoff))        
	    #else:
		plt.show()
    else:
	if 1:
	    path = "/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/result/ChapmanKolmogorovTest"
	    plt.figure()
	    plt.plot(range(len(Prob_tau)),Prob_tau)
	    plt.hold(True)
	    plt.plot(range(len(Prob_ktau)),Prob_ktau)
	    plt.legend(['Prob_tau','Prob_ktau'])
	    plt.title("Cutoff = %0.1f tau = %d "%(cutoff,tau))
	    plt.xlabel('K')
	    plt.ylabel('State Population')
	    plt.ylim(0,1)
	    
	    if 1:
                figname = "CKtest_tau%d_RMSDCutoff%0.1f.png"%(tau,cutoff)
                print "Save to %s"%(figname)
		plt.savefig(figname)        
	    #else:
		#plt.show()	    
	    
    
if __name__ == '__main__' :
    pass
