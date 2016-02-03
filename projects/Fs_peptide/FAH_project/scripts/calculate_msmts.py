import os,sys
import numpy as np
import mdtraj as md
from sklearn.externals import joblib
from itertools import combinations
from msmbuilder.featurizer import AtomPairsFeaturizer,RMSDFeaturizer
from msmbuilder.decomposition import tICA
from msmbuilder.cluster import KMeans,KCenters,KMedoids
from msmbuilder.msm import MarkovStateModel


verbose = True

ids = ['combined']+range(6383,6391)
mapid = 40

for pid in ids: 
    if pid == 'combined':
        a = np.load('Assignments.fixed.Map{}.npy'.format(mapid))
    else:
        a = np.load('Assignments-{}.fixed.Map{}.npy'.format(pid,mapid))
    if verbose:
        print "Building MSMs:"
    lagtimes = [1,10,20,30,40,50,100,150,200]
    msmts = []
    for lagtime in lagtimes:
        if verbose:
            print "\tLagtime: %d"%lagtime
        msm = MarkovStateModel(lag_time=lagtime).fit(a)
        msmts.append(msm.timescales_)
    output_dir = "Data-{}-macro{}".format(pid,mapid)
    lagtime_fn = os.path.join(output_dir,"lagtimes.txt")
    msmts_fn = os.path.join(output_dir,"ImpliedTimescales.npy")
    np.savetxt(lagtime_fn,lagtimes)
    np.save(msmts_fn,msmts)
    if verbose:
        print "Wrote: %s and %s"%(lagtime_fn,msmts_fn)

