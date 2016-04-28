import os,sys
import numpy as np
from msmbuilder import tpt
from scipy.sparse import coo_matrix
from scipy.io import mmwrite
from sklearn.externals import joblib

def rate_estimation_tpt(source_states,sink_states,msm):
    tau = 5
    #timescales_macro = 3788.1 #1st eigenmode
    #timescales_micro = 12095.2 #1st eigenmode, may not be correct!

    c_p53 = 7.1*10**-3 #Mol

    committors = tpt.committors(sources=source_states,sinks=sink_states,msm=msm)
    back_com = 1-committors

    #print committors
    F = 0
    for i in source_states:
        for j in range(msm.n_states_):
            if j not in source_states:
                F += msm.populations_[i]*msm.transmat_[i,j]*committors[j]

    kab = F/((tau*1e-9)*np.dot(msm.populations_,back_com))
    #print kab*(timescales_macro/timescales_micro)
    #print kab*(1/c_p53)

    #print "kAB:",kab*(1/c_p53)*(timescales_macro/timescales_micro)
    return kab

def rate_estimation_msm(assignments,mapping):
    

    return None


rmsd = np.load("rmsd-states-micro1000.npy")
d = np.load("mdm2-p53-distance-states-micro1000.npy")

rmsd_cutoff = 0.25 #nm
d_cutoff = 2.15 #nm

c_p53 = 7.1*10**-3 #M

d_native = d[860]
print "native distance:",d_native
delta = 0.05

fold_ndx = np.where(rmsd<rmsd_cutoff)[0]
unfold_ndx = np.where(rmsd>=rmsd_cutoff)[0]

if 0:
    bound_ndx = np.where(d < d_cutoff)[0]
    unbound_ndx = np.where(d >= d_cutoff)[0]
else:
    bound_ndx = np.where(abs(d-d_native)<delta)[0]
    unbound_ndx = np.where(abs(d-d_native)>=delta)[0]

region_dict = {"fold":fold_ndx,"unfold":unfold_ndx,
               "bound":bound_ndx,"unbound":unbound_ndx}


unfold_unbound = np.intersect1d(region_dict["unfold"],region_dict["unbound"])
fold_unbound = np.intersect1d(region_dict["fold"],region_dict["unbound"])
unfold_bound = np.intersect1d(region_dict["unfold"],region_dict["bound"])
fold_bound = np.intersect1d(region_dict["fold"],region_dict["bound"])


regions = {"unfold_unbound":unfold_unbound,"fold_unbound":fold_unbound,
           "unfold_bound":unfold_bound,"fold_bound":fold_bound}

msm = joblib.load("../MSMs-1000-microstates/MSMs-1000-microstates.pkl")

unfold_region = np.union1d(unbound_ndx,unfold_bound)

k_off_overall = rate_estimation_tpt(fold_bound,unfold_region,msm)
k_on_overall = rate_estimation_tpt(unfold_region,fold_bound,msm)
print "Overall k_off: %.3e 1/s"%k_off_overall
print "Overall k_on: %.3e 1/(M*s)"%(k_on_overall/(msm.populations_[fold_bound].sum()*c_p53))

for region in regions.keys():
    print region,"population:",msm.populations_[regions[region]].sum(),msm.populations_[regions[region]].sum()*c_p53,"mol/L"

for source_region in regions.keys():
    for sink_region in regions.keys():
        if sink_region == source_region:
            continue
        sink_states = regions[sink_region]
        source_states = regions[source_region]
        rate = rate_estimation_tpt(source_states,sink_states,msm)
        if "unbound" in source_region:
            rate1 = rate/(msm.populations_[source_states].sum()*c_p53)
            print "From %s to %s, transition rate=%.3e,rate/concentration = %.3e 1/(M*s)"%(source_region,sink_region,rate,rate1)
        print "From %s to %s, transition rate=%.3e"%(source_region,sink_region,rate)
