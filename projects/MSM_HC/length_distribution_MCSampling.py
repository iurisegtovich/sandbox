import os,sys
import numpy as np
sys.path.append('/Users/tud51931/voelzlab/analysis/LifsonRoig/src')
import LifsonRoigTools as LRTools
import random
import copy
import matplotlib.pyplot as plt






def random_change(sequence):
    sequence_new = copy.copy(sequence)
    i = random.randrange(len(sequence))
    sequence_new[i] = (sequence[i]+1)%2
    return sequence_new

steps = 100000

l = np.loadtxt('/Users/tud51931/voelzlab/analysis/LifsonRoig/scripts/test_Fs_RRR_ff03/Loglikelihood.dat')
w = np.loadtxt('/Users/tud51931/voelzlab/analysis/LifsonRoig/scripts/test_Fs_RRR_ff03/w_params.dat')
v = np.loadtxt('/Users/tud51931/voelzlab/analysis/LifsonRoig/scripts/test_Fs_RRR_ff03/v_params.dat')

I = np.argsort(l)
w_max = w[I[-1]]
v_max = v[I[-1]]

sequence = np.zeros(len(w_max))
a = LRTools.LifsonRoigTools()
i = 0
count = np.zeros(len(w_max)+1)
trajectory = []
trajectory.append(sequence)
count[0] = count[0] + 1
accept = 0
while i <=steps:
    a.load_from_HCarray(sequence)
    a.calculate_w_v_array()
    weight = a.calculate_weight_frame(a.w_array,a.v_array,w_max,v_max)

    sequence_new = random_change(sequence)
    a.load_from_HCarray(sequence_new)
    a.calculate_w_v_array()
    weight_new = a.calculate_weight_frame(a.w_array,a.v_array,w_max,v_max)
    N_helix = a.w_array.sum()+a.v_array.sum()
    if random.random() <= min(1,weight_new/weight):
        print "Accept!-->P = %f"%(weight_new/weight)
        #print "Old Array:",sequence
        #print "New Array:",sequence_new
        accept = accept + 1.0
        sequence = sequence_new
        count[N_helix] = count[N_helix] + 1
        trajectory.append(sequence)
    i = i + 1
print "Done! Accept ratio:%f"%(accept/steps)
print count
plt.figure()
plt.plot(range(len(count)),count/np.sum(count),'o')
plt.show()


