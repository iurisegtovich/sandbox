import os,sys
import numpy as np

def rebuild_tProb(n_marcro,tProb, mapping):
    tProb_new = np.zeros((n_macro,n_macro))
    if isinstance(mapping,dict):
        keys,values = mapping.keys(),mapping.values()
        keys = np.array(keys)
        values = np.array(values)
        index_grid_keys = np.ix_(keys,keys)
        index_grid_values = np.ix_(values,values)
        tProb_new[index_grid_keys] = tProb[index_grid_values]
    return tProb_new
        

def run(n_macro=40, outdir = "tProb_rebuild"):
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    pid_list = range(6383,6391)

    msm_fn = ""
    combine_mBACE_fn = "Output_BACE+/combination%d.dat"%n_macro
    combining_mBACE = np.loadtxt(combine_mBACE_fn,dtype=int)
    mapping_mBACE_fn = "Output_BACE+/map%d.dat"%n_macro
    mapping_mBACE = np.loadtxt(mapping_mBACE_fn,dtype=int)

    print "uncombined:", mapping_mBACE[combining_mBACE==1]
    print "combined:", mapping_mBACE[combining_mBACE==0]

    tCounts_combined = np.zeros((n_macro,n_macro))
    tCounts_all_proj = []
    mapping_all_proj = []

    for pid in pid_list:
        assignments_fn = "../results/Assignments_BACE+/Assignments-{}.fixed.Map{}.BACE+.npy".format(pid,n_macro)
        assignments = np.load(assignments_fn)
        tCounts, mapping = core._transition_counts(assignments, lag_time=50, sliding_window=True)
        if len(mapping) < n_macro:
            tCounts_new = np.zeros((n_macro,n_macro))
            keys,values = mapping.keys(),mapping.values()
            keys = np.array(keys)
            values = np.array(values)
            #print keys
            #print values
            tCounts_new[keys.reshape(-1,1),keys] = tCounts[values.reshape(-1,1),values]
            tCounts_all_proj.append(tCounts_new)
            print mapping
        else:
            tCounts_all_proj.append(tCounts)
        mapping_all_proj.append(mapping)
        print pid,len(mapping)

