import os,sys
import numpy as np
import pickle
from sklearn.externals import joblib
from scipy.sparse import coo_matrix
from scipy.io import mmwrite
from msmbuilder.msm import MarkovStateModel


verbose = True

map_id = 40

for p_id in range(6383,6391): 
    assignments = np.load('Assignments-%d.fixed.Map%d.npy'%(p_id,map_id))

    if verbose:
        print "Building MSMs:"
    lagtime = 50
    msm = MarkovStateModel(lag_time=lagtime).fit(assignments)
    print msm
    model_dir = "MSMs-%d-macro%d"%(p_id,map_id)
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)
    msm_fn = os.path.join(model_dir,"MSMs-%d-macro%d.pkl"%(p_id,map_id))
    if os.path.exists(msm_fn):
        print "%s already exists! Quit..."%msm_fn
        sys.exit()
    joblib.dump(msm,msm_fn)
    print "Wrote:%s"%msm_fn

    output_dir = "Data-%d-macro%d"%(p_id,map_id)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    fn_countsmat = os.path.join(output_dir,"tCounts.npy")
    fn_transmat = os.path.join(output_dir,"tProb.npy")
    fn_populations = os.path.join(output_dir,"Populations.dat")
    fn_mapping = os.path.join(output_dir,"Mapping.dat")
    fn_countsmat_sparse = os.path.join(output_dir,"tCounts.mtx")
    fn_transmat_sparse = os.path.join(output_dir,"tProb.mtx")

    np.save(fn_countsmat,msm.countsmat_)
    np.save(fn_transmat,msm.transmat_)
    np.savetxt(fn_populations,msm.populations_)
    mapping = []
    for i in range(len(msm.mapping_)):
        try:
            mapping.append(msm.mapping_[i])
        except KeyError:
            mapping.append(-1)
    np.savetxt(fn_mapping,mapping,fmt="%d")

    mmwrite(fn_countsmat_sparse,coo_matrix(msm.countsmat_))
    mmwrite(fn_transmat_sparse,coo_matrix(msm.transmat_))


