import os, sys
import numpy as np
from msmbuilder.msm import MarkovStateModel
import matplotlib.pyplot as plt


def CalculateStatesProbability(T_tau, T_ktau, populations_tau, populations_ktau, mapping_tau, mapping_ktau,
                               k, states):
    probability_tau = 0
    probability_ktau = 0
    populations_tau_i = np.zeros(len(populations_tau))
    populations_ktau_i = np.zeros(len(populations_ktau))
    for i in states:
        try:
            if mapping_tau.has_key(i) and mapping_ktau.has_key(i):
                populations_tau_i[i] = populations_tau[int(mapping_tau[i])]
                populations_ktau_i[i] = populations_ktau[int(mapping_ktau[i])]
            else:
                print "Warning:states %d does not match, remove state %d in the list" % (i, i)
                states.remove(i)
                if len(states) == 0:
                    return -1,-1
        except IndexError:
            print "state number %d is out of bound, remove state %d in the list" % (i, i)
            states.remove(i)
            if len(states) == 0:
                return -1,-1
    # Normalize the initial populations
    populations_tau_i /= populations_tau_i.sum(0)
    populations_ktau_i /= populations_ktau_i.sum(0)
    p_tau = np.matrix(populations_tau_i)
    for _ in range(k):
        p_tau = p_tau * T_tau
    p_ktau = np.matrix(populations_ktau_i) * T_ktau
    for i in states:
        probability_tau += p_tau[0, int(mapping_tau[i])]
        probability_ktau += p_ktau[0, int(mapping_ktau[i])]

    return probability_tau, probability_ktau


def ChapmanKolmogorovTest(assignments, klist = [1,2,3,4,5], lagtime=50, states=None):
    msm = MarkovStateModel(lag_time=lagtime, n_timescales=10)
    msm.fit(assignments)
    p_tau = msm.populations_
    T_tau = msm.transmat_
    mapping_tau = msm.mapping_

    prob_tau_all = []
    prob_ktau_all = []


    if states == "all" or states is None:
        states = range(len(p_tau))

    for k in klist:
        lagtime_long = k * lagtime
        print "long lagtime:",lagtime_long
        msm = MarkovStateModel(lag_time=lagtime_long, n_timescales=10)
        msm.fit(assignments)
        p_ktau = msm.populations_
        T_ktau = msm.transmat_
        mapping_ktau = msm.mapping_
        probability_tau, probability_ktau = CalculateStatesProbability(T_tau, T_ktau, p_tau, p_ktau,
                                                                       mapping_tau, mapping_ktau, k, states)

        prob_tau_all.append(probability_tau)
        prob_ktau_all.append(probability_ktau)

    return prob_tau_all, prob_ktau_all


def test(assignments_fn):
    assignments = np.load(assignments_fn)
    prob_tau_all, probl_ktau_all = ChapmanKolmogorovTest(assignments, klist=range(1,101), lagtime=50, states = [20])
    print prob_tau_all
    print probl_ktau_all

def plot_CK_test(assignments_fn, sequence = "combined",stateid=20):
    klist=range(1,101)
    #klist.remove(5)
    #klist.remove(49)
    assignments = np.load(assignments_fn)
    prob_tau_all, prob_ktau_all = ChapmanKolmogorovTest(assignments, klist=klist, lagtime=50, states = [stateid])
    prob_tau_all = [1.] + prob_tau_all
    prob_ktau_all = [1.] + prob_ktau_all
    prob_tau_all=np.array(prob_tau_all)
    prob_ktau_all=np.array(prob_ktau_all)
    lagtime_list = np.array([0]+klist)*5
    print len(prob_tau_all),len(lagtime_list)
    posInd = np.where(prob_tau_all>=0)[0]
    plt.figure()
    plt.title("Fs-%s State %d"%(sequence,stateid))
    plt.ylabel("State residence probability")
    plt.ylim([0,1])
    plt.plot(lagtime_list[posInd],prob_tau_all[posInd],'r--')
    plt.plot(lagtime_list[posInd],prob_ktau_all[posInd],'b')
    figure_fn = "CKtest-Fs-%s-State%d.png"%(sequence,stateid)
    plt.savefig(figure_fn)
    print "Saved: %s"%figure_fn

def plot_all():
    sequences = ["EEE","EER","ERE","ERR","REE","RER","RRE","RRR"]
    for i,projid in enumerate(range(6383,6391)):
        seqid = sequences[i]
        assignments_fn = "../test/Assignments-%d.fixed.Map40.npy"%projid
        #test(assignments_fn)
        for stateid in [13,20,11,2]:
            plot_CK_test(assignments_fn, sequence = seqid, stateid=stateid)

def plot_CK_test_RRR(assignments_fn, sequence = "combined",stateid=20):
    klist=range(1,101)
    klist.remove(5)
    klist.remove(49)
    assignments = np.load(assignments_fn)
    prob_tau_all, prob_ktau_all = ChapmanKolmogorovTest(assignments, klist=klist, lagtime=50, states = [stateid])
    prob_tau_all = [1.] + prob_tau_all
    prob_ktau_all = [1.] + prob_ktau_all
    prob_tau_all=np.array(prob_tau_all)
    prob_ktau_all=np.array(prob_ktau_all)
    lagtime_list = np.array([0]+klist)*5
    print len(prob_tau_all),len(lagtime_list)
    posInd = np.where(prob_tau_all>=0)[0]
    plt.figure()
    plt.title("Fs-%s State %d"%(sequence,stateid))
    plt.ylabel("State residence probability")
    plt.ylim([0,1])
    plt.plot(lagtime_list[posInd],prob_tau_all[posInd],'r--')
    plt.plot(lagtime_list[posInd],prob_ktau_all[posInd],'b')
    figure_fn = "CKtest-Fs-%s-State%d.png"%(sequence,stateid)
    plt.savefig(figure_fn)
    print "Saved: %s"%figure_fn


if __name__ == "__main__":
     assignments_fn = "../test/Assignments-6390.fixed.Map40.npy"
     plot_CK_test_RRR(assignments_fn, sequence = "RRR", stateid=20)
     #plot_all()
