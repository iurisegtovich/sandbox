# GMRQ hyperparameter selection
import os,sys
import pandas as pd
import numpy as np
import pickle
from msmbuilder.cluster import KMeans
from sklearn.cross_validation import KFold
from msmbuilder.msm import MarkovStateModel

def sub_sampling_data(raw_data,stride):
    reduced_data = []
    for traj in raw_data:
        reduced_data.append(traj[::stride])
    return reduced_data

if len(sys.argv) < 2:
    print "Usage: python this.py n_tica_component"
    sys.exit()
n_components = [int(sys.argv[1])]
#n_components=[2,4,6,8,10,15]
nFolds=5
for n_tIC in n_components:
    tica_fn = "tICA_data_lagtime50_tIC%d.npy"%n_tIC
    if not os.path.exists(tica_fn):
        print "%s not exists!"%tica_fn
        continue

    tica_data = np.load(tica_fn)

    results = []

    n_clusters = [100,200,400,600,800,1000,1200,1500,2000,2500,3000]

    #n_clusters = [1200,1500,2000]

    #n_clusters = [3500,3500,4000,4500,5000,6000]
    lagtime = 50

    for n in n_clusters:
        kmeans = KMeans(n_clusters=n,n_jobs=-1)
        print "Clustering data to %d clusters..."%n
        for fold in range(nFolds):
            train_data=[]
            test_data=[]
            for i in range(len(tica_data)):
                cv = KFold(len(tica_data[i]),n_folds=nFolds)
        	for current_fold,(train_index,test_index) in enumerate(cv):
                    if current_fold == fold:
              	        train_data.append(tica_data[i][train_index])
        	        test_data.append(tica_data[i][test_index])
            reduced_train_data = sub_sampling_data(train_data,stride=100)
            kmeans.fit(reduced_train_data)
            assignments_train = kmeans.predict(train_data)
            assignments_test = kmeans.predict(test_data)
            msm = MarkovStateModel(lag_time=lagtime)
            msm.fit(assignments_train)
            train_score = msm.score_
            test_score = msm.score(assignments_test)

            results.append({'train_score':train_score,
                            'test_score':test_score,
                            'n_states':n,
                            'fold':fold,
                            'timescales':msm.timescales_})

        results = pd.DataFrame(results)
        print results
        output_fn = "GMRQ_MSMs_score_for_tica_n_%d_%d.pkl"%(n_tIC,n)
        with open(output_fn,'wb') as result_fn:
            pickle.dump(results,result_fn)
            result_fn.close()
