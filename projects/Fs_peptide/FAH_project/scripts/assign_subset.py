import os,sys
import numpy as np
import mdtraj as md
from sklearn.externals import joblib


verbose = True

tica_data=np.load('../../ticas_n_8.npy')
with open('../../trajectories.log') as log:
    trajs_log = log.readlines()
kmeans = joblib.load('kmeans_model_n_1200/kmeans-combined.pkl')
print len(tica_data)

for p_id in range(6383,6391):
    p_index = [i for i,j in enumerate(trajs_log) if str(p_id) in j]
    print p_index
    assignments = kmeans.predict(tica_data[p_index])
    assignments_fn = "Assignments-%d.npy"%p_id
    np.save(assignments_fn,assignments)
    if verbose:
        print "Wrote: %s"%assignments_fn
