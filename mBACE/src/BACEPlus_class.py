#BACE Plus
import os
import sys
import numpy as np
import argparse
import functools
import multiprocessing
from scipy.io import mmread
import scipy.sparse

def check_matrices_shape(cMats):

    s = cMats[0].shape
    for i in range(len(cMats)):
        if cMats[i].shape != s:
            raise TypeError("Matrix %d has a different shape!"%i)

class BacePlus():

    def __init__(self, cMats, n_macro=2, nProc=1):

        check_matrices_shape(cMats)

        self.n_models = len(cMats)
        self.n_macro = 2
        self.n_micro = cMats[0].shape[0]

        self.uncombined = np.ones(self.n_micro, dtype=np.int32)
        self.stateskeep = np.arange(self.n_micro, dtype=np.int32)
        self.unmerged = np.ones(self.n_micro, dtype=np.int32)
        self.indices = np.arange(self.n_micro, dtype=np.int32)
        self.map = np.arage(self.n_micro, dtype=np.int32)
        self.pseud = np.ones(self.n_micro, dtype=np.float32)/self.n_micro
        self.nProc = nProc


    def preCombineSparse_BacePlus(self):

        n_micro = cMats[0].shape[0]
        w_multi = np.zeros(((cMats),self.n_micro))
        for i,c in enumerate(cMats):
            w_multi[i][:] = np.array(c.sum(axis=1)).flatten()
        w_multi += 1

        map = np.arange(n_micro,dtype=np.int32)
        pseud = np.ones(n_micro,dtype=np.float32)/n_micro

        indices = np.arange(n_micro, dtype=np.int32)
        statesKeep = np.arange(n_micro, dtype=np.int32)
        statesCombine = np.arange(n_micro, dtype=np.int32)
        uncombined = np.ones(n_micro, dtype=np.int8)

        nInd = len(indices)
        # Compute the JSD in the same matrix
        if nInd >1 and nProc >1:
            if nInd < nProc:
                nProc=nInd
            pool = multiprocessing.Pool(processes=nProc)
            stepSize = int(nInd/nProc)
            if nInd % stepSize > 3:
                dlims = zip(range(0, nInd, stepSize),
                range(stepSize, nInd, stepSize) + [nInd])
            else:
                dlims = zip(range(0, nInd - stepSize, stepSize),
                range(stepSize, nInd - stepSize, stepSize) + [nInd])

            args = []
            for start, stop in dlims:
                args.append(indices[start:stop])
            print args
            result = pool.map_async(
                functools.partial(surprisalSparseHelper, cMats=cMats, w_multi=w_multi, statesKeep=statesKeep, uncombined=uncombined), args)
            result.wait()
            d = np.concatenate(result.get())
            pool.close()
        else:
            d = surprisalSparseHelper(
                indices, cMats, w_multi, statesKeep, uncombined)

        return d

def calSurprisalSparse(indices, cMats, w_multi, stateKeep, uncombined):
    d = np.zeros(indices.shape[0], dtype=np.float32)
    for i in range(indices.shape[0]):
        cRows = [c[i] for c in range(len(cMats))]


def surprisalSparseHelper(indices, cMats, w_multi, statesKeep, uncombined):
    d = np.zeros(indices.shape[0], dtype=np.float32)
    for i,index in enumerate(indices):
        cRows = np.array([c[index].toarray()[0]+ 1.0 /c.shape[0] for c in cMats])
        wRows = np.array([w[index] for w in w_multi])
        #print index,cRows
        p_tot = cRows.sum(axis=0)/wRows.sum(axis=0)
        for j in range(cRows.shape[0]):
            p_row = cRows[j]/wRows[j]
            d[i] += cRows[j].dot(np.log(p_row/p_tot))
    return d

def test_precombine():
    cMats = []
    for i in range(6383,6391):
        #fn = "/Users/gfzhou/git/sandbox/bace+/test/tCounts/Fs-%d-tCounts-macro40.mtx"%i
        fn = "tCounts-micro-%d.mtx"%i
        cMats.append(mmread(fn).tolil())
    a = BacePlus(cMats=cMats, nProc=4)
    d = preCombineSparse_BacePlus(cMats,nProc=4)
    print d

if __name__ == "__main__":
    test_precombine()




