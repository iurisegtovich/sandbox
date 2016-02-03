import os,sys
import numpy as np
from sklearn.externals import joblib
import matplotlib.pyplot as plt

map_id = 40
pids = range(6383,6391)
mut = ['EEE','EER','ERE','ERR','REE','RER','RRE','RRR']
p = 1./40.*np.ones(40)

for i,p_id in enumerate(pids):
    
    model_dir = "MSMs-%d-macro%d"%(p_id,map_id)
    if not os.path.exists(model_dir):
        print "%s doesn't exist! Exit!"%model_dir
    msm_fn = os.path.join(model_dir,"MSMs-%d-macro%d.pkl"%(p_id,map_id))
    msm = joblib.load(msm_fn)
    for j in range(5):
        print mut[i],j,np.dot(p,msm.right_eigenvectors_[:,j])





