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

tica_data=np.load('../../ticas_n_8.npy')

reduced_data = []
for i in range(len(tica_data)):
    reduced_data.append(tica_data[i][::100,:])

if verbose:
    print "Clustering."
kmeans = KMeans(n_clusters=1200).fit(reduced_data)
Gen_fn = "Gens.npy"
np.save(Gen_fn,kmeans.cluster_centers_)
if verbose:
    print "Wrote: %s"%Gen_fn
model_dir = "kmeans_model_n_1200"
if not os.path.exists(model_dir):
    os.makedirs(model_dir)
model_fn = os.path.join(model_dir,'kmeans-combined.pkl')
joblib.dump(kmeans,model_fn)
if verbose:
    print "Saved cluster model to %s"%model_fn
if verbose:
    print "Assigning.."
assignments = kmeans.predict(tica_data)
if verbose:
    print "Wrote assignments"
np.save('Assignments.npy',assignments)

if verbose:
    print "Building MSMs:"
lagtimes = [1,10,20,30,40,50,100,150,200]
msmts = []
for lagtime in lagtimes:
    if verbose:
        print "\tLagtime: %d"%lagtime
    msm = MarkovStateModel(lag_time=lagtime).fit(assignments)
    msmts.append(msm.timescales_)
lagtime_fn = "lagtimes.txt"
msmts_fn = "ImpliedTimescales.npy"
np.savetxt(lagtime_fn,lagtimes)
np.save(msmts_fn,msmts)
if verbose:
    print "Wrote: %s and %s"%(lagtime_fn,msmts_fn)

