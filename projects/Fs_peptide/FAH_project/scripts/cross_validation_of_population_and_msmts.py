import os,sys
import numpy as np
import pickle
import pandas as pd
from sklearn.cross_validation import KFold
from sklearn.externals import joblib
from scipy.sparse import coo_matrix
from scipy.io import mmwrite
from msmbuilder.msm import MarkovStateModel


def case1():
    map_id = 40
    for p_id in range(6383,6391):

        assignments = np.load('Assignments-%d.fixed.Map%d.npy'%(p_id,map_id))
        cv=KFold(len(assignments),n_folds=10)
        lagtime = 50
        msm = MarkovStateModel(lag_time=lagtime)
        pops = []
        msmts = []
        for fold,(train_index,test_index) in enumerate(cv):
            assignments_train = assignments[train_index]
            msm.fit(assignments_train)
            if len(msm.populations_)==40:
                pops.append(msm.populations_)

            msmts.append(msm.timescales_)


        output_dir = "Data-%d-macro%d"%(p_id,map_id)
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        fn_populations = os.path.join(output_dir,"Populations-10fold.npy")
        fn_msmts = os.path.join(output_dir,"ImpliedTimescales-10fold.npy")

        np.save(fn_populations,pops)
        np.save(fn_msmts,msmts)
        print "Saved: {},{}".format(fn_populations,fn_msmts)

def case2_micro_combined():
    assignments = np.load('Assignments.npy')
    lagtime = 50
    msmts = []
    msm = MarkovStateModel(lag_time=lagtime)
    cv=KFold(len(assignments),n_folds=10)
    for fold,(train_index,test_index) in enumerate(cv):
        assignments_train = assignments[train_index]
        msm.fit(assignments_train)
        msmts.append(msm.timescales_)
    fn_msmts = os.path.join('Data-combined',"ImpliedTimescales-10fold.npy")
    np.save(fn_msmts,msmts)
    print "Saved: {}".format(fn_msmts)

def case3_macro_combined():
    assignments = np.load('Assignments.fixed.Map40.npy')
    lagtime = 50
    msmts = []
    msm = MarkovStateModel(lag_time=lagtime)
    cv=KFold(len(assignments),n_folds=10)
    for fold,(train_index,test_index) in enumerate(cv):
        assignments_train = assignments[train_index]
        msm.fit(assignments_train)
        msmts.append(msm.timescales_)
    fn_msmts = os.path.join('Data-combined-macro40',"ImpliedTimescales-10fold.npy")
    np.save(fn_msmts,msmts)
    print "Saved: {}".format(fn_msmts)


if __name__ == "__main__" :
    case2_micro_combined()
    case3_macro_combined()

