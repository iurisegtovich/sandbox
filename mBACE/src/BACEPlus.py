#BACE Plus
import os
import sys
import numpy as np
from copy import copy
import argparse
import functools
import multiprocessing
from scipy.io import mmread
import scipy.sparse

def check_cMats(cMats):

    s = cMats[0].shape
    for i in range(len(cMats)):
        if cMats[i].shape != s:
            raise TypeError("Matrix %d has a different shape!"%i)

def renumberMap(map, stateDrop):
    for i in xrange(map.shape[0]):
        if map[i] >= stateDrop:
            map[i] -= 1
    return map

def getIndices(cMats, combinedStateInds, chunkSize, updateSingleState=None):

    if scipy.sparse.issparse(cMats[0]):
        cTotal = scipy.sparse.lil_matrix(cMats[0].shape)
    else:
        cTotal = np.zeros(cMats[0].shape)
    for i in range(len(cMats)):
        cTotal += cMats[i]
    isSparse = scipy.sparse.issparse(cTotal)

    indices = []

    stateInds = combinedStateInds

    if updateSingleState != None:
        if isinstance(updateSingleState,int):
            stateInds = [updateSingleState]
        elif isinstance(updateSingleState,list):
            stateInds = updateSingleState
    print "stateInds",stateInds
    for s in stateInds:

        if isSparse:
            positive = set(np.where(cTotal[s,:].toarray()[0]>0)[0])
        else:
            positive = set(np.where(cTotal[s,:]>0)[0])

        if updateSingleState != None:
            combinedInds = set(combinedStateInds[np.where(combinedStateInds!= updateSingleState)[0]])
            dest = np.array(list(positive&combinedInds))
            dest = np.array(list(dest))
        else:
            combinedInds = set(combinedStateInds[np.where(combinedStateInds > s)[0]])
            dest = np.array(list(positive&combinedInds))
            dest = np.array(list(dest))

        if dest.shape[0] == 0:
            continue
        elif dest.shape[0] < chunkSize:
            indices.append((s, dest))
        else:
            i = 0
            while dest.shape[0] > i:
                if i + chunkSize > dest.shape[0]:
                    indices.append((s, dest[i:]))
                else:
                    indices.append((s, dest[i:i + chunkSize]))
                i += chunkSize
    return indices

def splitIndices(indices, nProc):

    nInds = len(indices)
    if nInds < nProc:
        nProc = nInds

    stepSize = int(nInds / nProc)
    if nInds % stepSize > 3:
        dlims = zip(
            range(0, nInds, stepSize), range(stepSize, nInds, stepSize) + [nInds])
    else:
        dlims = zip(range(0, nInds - stepSize, stepSize),
                    range(stepSize, nInds - stepSize, stepSize) + [nInds])
    args = []
    for start, stop in dlims:
        args.append(indices[start:stop])

    return args


def run(cMats, nMacro, nProc, preCombine, multiDist, calcSurprisal, chunkSize=100):

    check_cMats(cMats)
    nModels = len(cMats)
    nMicro = cMats[0].shape[0]
    nCurrentStates = cMats[0].shape[0]
    w_multi = np.zeros((len(cMats),nMicro))
    for i,c in enumerate(cMats):
        w_multi[i][:] = np.array(c.sum(axis=1)).flatten()
    w_multi += 1
    wTotal = w_multi.sum(axis=0)

    statesCombined, uncombined= preCombine(cMats,nProc)
    print statesCombined

    map = np.arange(nMicro, dtype=np.int32)
    indices = np.arange(nMicro, dtype=np.int32)
    statesKeep = np.arange(nMicro, dtype=np.int32)
    unmerged = np.ones((nModels,nMicro), dtype=np.int8)

    newInd = getIndices(cMats, statesCombined, chunkSize)

    print "newInd",newInd
    if scipy.sparse.issparse(cMats[0]):
        dMat = scipy.sparse.lil_matrix(cMats[0].shape)
    else:
        dMat = np.zeros(cMats[0].shape, dtype=np.float32)

    sMat = np.zeros(nMicro, dtype=np.float32)

    outDir = '.'
    fBayesFact = open("%s/bayesFactors.dat" % outDir, 'w')
    dMat, minX1, minY1, maxD = calcCombinedDMat(
        cMats, w_multi, fBayesFact, newInd, dMat, nProc, statesKeep, multiDist, unmerged[0], chunkSize)
    #print "dMat,minX1,minY1,maxD", dMat, minX1, minY1, maxD

    sMat, minX2, maxS = calcSMat(
        cMats, w_multi, fBayesFact, sMat, nProc, statesKeep, calcSurprisal, uncombined, chunkSize)
    #print "sMat, minX2, maxS", sMat, minX2, maxS
    print "sMat",sMat

    i = 0

    while nCurrentStates > nMacro:
        print "maxD, maxS:", maxD, maxS

        if maxD > maxS:
            print "Iteration %d, lumping two closest combined states: %d,%d"%(i,minX1,minY1)
            print "nCurrentState:%d, nCombinedState:%d"%(nCurrentStates,len(np.where(uncombined!=1)[0]))
            cMats, w_multi, newInd, dMat, map, statesKeep, unmerged, minX1, minY1, maxD = lumpTwoStates(
                cMats, w_multi, fBayesFact, newInd, dMat, nProc, map, statesKeep, minX1, minY1, multiDist,
                unmerged, uncombined, chunkSize)
            nCurrentStates -= 1
        else:
            print "Iteration %d, combining states %d"%(i,minX2)
            print "nCurrentState:%d, nCombinedState:%d"%(nCurrentStates,len(np.where(uncombined!=1)[0]))
            sMat, dMat, minX1, minY1, minX2, maxS, maxD, uncombined = combineTwoStates(
                cMats, w_multi, fBayesFact, sMat, dMat, nProc, statesKeep, calcSurprisal,
                multiDist, minX2, uncombined, unmerged, chunkSize)
        i += 1


def combineTwoStates(cMats, w_multi, fBayesFact, sMat, dMat, nProc, statesKeep, calcSurprisal, multiDist, minX,
                     uncombined, unmerged, chunkSize, record=False):

    #doulbe check that the state hasn't been combined
    if not uncombined[minX]:
        raise ValueError("State %d already combined, check bugs...."%minX)
    uncombined[minX] = 0
    sMat, minX2, maxS = calcSMat(
        cMats, w_multi, fBayesFact, sMat, nProc, statesKeep, calcSurprisal, uncombined, chunkSize, updateSingleState=minX)

    statesCombined = np.where(uncombined==0)[0]

    indRecalc = getIndices(
        cMats, statesCombined, chunkSize, updateSingleState=int(minX))

    dMat, minX, minY, maxD = calcCombinedDMat(
        cMats, w_multi, fBayesFact, indRecalc, dMat, nProc, statesKeep, multiDist, unmerged[0], chunkSize, record)

    return sMat, dMat, minX, minY, minX2, maxS, maxD, uncombined






def lumpTwoStates(cMats, w_multi, fBayesFact, indRecalc, dMat, nProc, map, statesKeep, minX, minY, multiDist,
                  unmerged, uncombined, chunkSize, record=True):

    #double check before lumping two states
    if uncombined[minX] != 0 or uncombined[minY] != 0:
        raise ValueError("Lumping before combined, check bugs...")
    # Not sure whether to add psued counts or not...
    for i in range(len(cMats)):
        cIsSparse = scipy.sparse.issparse(cMats[i])
        if unmerged[i][minX]:
            cMats[i][minX, statesKeep] += unmerged[i][statesKeep] * 1.0 / cMats[i].shape[0]
            unmerged[i][minX] = 0
            if cIsSparse:
                cMats[i] = cMats[i].tolil()
                cMats[i][statesKeep, minX] += np.matrix(
                    unmerged[i][statesKeep]).transpose() * 1.0 / cMats[i].shape[0]
            else:
                cMats[i][statesKeep, minX] += unmerged[i][statesKeep] * 1.0 / cMats[i].shape[0]

        if unmerged[i][minY]:
            cMats[i][minY, statesKeep] += unmerged[i][statesKeep] * 1.0 / cMats[i].shape[0]
            unmerged[i][minY] = 0
            cIsSparse = scipy.sparse.issparse(cMats[i])
            if cIsSparse:
                cMats[i] = cMats[i].tolil()
                cMats[i][statesKeep, minY] += np.matrix(
                    unmerged[i][statesKeep]).transpose() * 1.0 / cMats[i].shape[0]
            else:
                cMats[i][statesKeep, minY] += unmerged[statesKeep] * 1.0 / cMats[i].shape[0]
        cMats[i][minX, statesKeep] += cMats[i][minY, statesKeep]
        cMats[i][statesKeep, minX] += cMats[i][statesKeep, minY]
        cMats[i][minY, statesKeep] = 0
        cMats[i][statesKeep, minY] = 0
        dMat[minX,:] = 0
        dMat[:, minX] = 0
        dMat[minY,:] = 0
        dMat[:, minY] = 0
        if cIsSparse:
            cMats[i] = cMats[i].tocsr()
        w_multi[i][minX] += w_multi[i][minY]
        w_multi[i][minY] = 0
    statesKeep = statesKeep[np.where(statesKeep != minY)[0]]
    print "minY",minY
    uncombined[minY] = -1
    indChange = np.where(map == map[minY])[0]
    map = renumberMap(map, map[minY])
    map[indChange] = map[minX]
    statesCombined = np.where(uncombined==0)[0]
    print "statesCombined:",statesCombined

    indRecalc = getIndices(
        cMats, statesCombined, chunkSize, updateSingleState=int(minX))
    dMat, minX, minY, maxD = calcCombinedDMat(
        cMats, w_multi, fBayesFact, indRecalc, dMat, nProc, statesKeep, multiDist, unmerged[0], chunkSize, record)

    return cMats, w_multi, indRecalc, dMat, map, statesKeep, unmerged, minX, minY, maxD


def calcSMat(cMats, w_multi, fBayesFact, sMat, nProc, statesKeep, calcSurprisal, uncombined, chunkSize, updateSingleState=None):

    if updateSingleState != None:
        if isinstance(updateSingleState,int):
            sMat[updateSingleState] = 0.0
            minX = sMat.argmax()
            maxS = sMat[minX]
            return sMat, minX, maxS
        else:
            raise TypeError("updateSingleState should be int type, get a %s instead"%(type(updateSingleState)))


    indUncombined = np.where(uncombined == 1)[0]

    nInd = len(indUncombined)
    if nInd > 1 and nProc > 1:
        pool = multiprocessing.Pool(processes=nProc)
        args = splitIndices(indUncombined, nProc)
        result = pool.map_async(
            functools.partial(calcSurprisal, cMats=cMats, w_multi=w_multi, statesKeep=statesKeep, uncombined=uncombined),args)
        result.wait()
        d = np.concatenate(result.get())
        pool.close()
    else:
        d = calcSurprisal(indUncombined, cMats, w_multi, statesKeep, uncombined)
    for i in xrange(len(indUncombined)):
        sMat[indUncombined[i]] = d[i]
    minX = sMat.argmax()
    maxS = sMat[minX]

    return sMat, minX, maxS



def calcCombinedDMat(cMats, w_multi, fBayesFact, indRecalc, dMat, nProc, statesKeep, multiDist, unmerged, chunkSize,
                     record=True):

    if scipy.sparse.issparse(cMats[0]):
        cTotal = scipy.sparse.lil_matrix(cMats[0].shape)
    else:
        cTotal = np.zeros(cMats[0].shape)
    for i in range(len(cMats)):
        cTotal += cMats[i]
    wTotal = w_multi.sum(axis=0)

    nRecalc = len(indRecalc)
    if nRecalc > 1 and nProc > 1:
        pool = multiprocessing.Pool(processes=nProc)
        args = splitIndices(indRecalc, nProc)
        result = pool.map_async(
            functools.partial(
                multiDist, c=cTotal, w=wTotal, statesKeep=statesKeep, unmerged=unmerged, chunkSize=chunkSize), args)
        result.wait()
        d = np.vstack(result.get())
        pool.close()
    else:
        d = multiDist(indRecalc, c=cTotal, w=wTotal, statesKeep=statesKeep, unmerged=unmerged, chunkSize=chunkSize)
    for i in xrange(len(indRecalc)):
        dMat[indRecalc[i][0], indRecalc[i][1]] = d[i][:len(indRecalc[i][1])]

    #print "dMat",dMat

    # BACE BF inverted so can use sparse matrices
    if scipy.sparse.issparse(dMat):
        minX = minY = -1
        maxD = 0
        for x in statesKeep:
            if len(dMat.data[x]) == 0:
                continue
            pos = np.argmax(dMat.data[x])
            if dMat.data[x][pos] > maxD:
                maxD = dMat.data[x][pos]
                minX = x
                minY = dMat.rows[x][pos]
    else:
        indMin = dMat.argmax()
        maxD = dMat[indMin]
        minX = np.floor(indMin / dMat.shape[1])
        minY = indMin % dMat.shape[1]
    if record:
        fBayesFact.write("%d %f\n" %
                         (statesKeep.shape[0] - 1, 1. / dMat[minX, minY]))
    return dMat, minX, minY, maxD


def preCombineSparse(cMats, nProc):

    n_micro = cMats[0].shape[0]
    w_multi = np.zeros((len(cMats),n_micro))
    for i,c in enumerate(cMats):
        w_multi[i][:] = np.array(c.sum(axis=1)).flatten()
    w_multi += 1

    indices = np.arange(n_micro, dtype=np.int32)
    statesKeep = np.arange(n_micro, dtype=np.int32)
    uncombined = np.ones(n_micro, dtype=np.int8)

    nInd = len(indices)

    if nInd >1 and nProc >1:
        pool = multiprocessing.Pool(processes=nProc)
        args = splitIndices(indices, nProc)
        result = pool.map_async(
            functools.partial(surprisalSparseHelper, cMats=cMats, w_multi=w_multi, statesKeep=statesKeep, uncombined=uncombined), args)
        result.wait()
        d = np.concatenate(result.get())
        pool.close()
    else:
        d = surprisalSparseHelper(
            indices, cMats, w_multi, statesKeep, uncombined)

    statesCombined = np.where(d<1.1)[0]
    uncombined[statesCombined] = 0

    return statesCombined, uncombined

def multiDistSparse(indicesList, c, w, statesKeep, unmerged, chunkSize):
    d = np.zeros((len(indicesList), chunkSize), dtype=np.float32)
    for j in xrange(len(indicesList)):
        indices = indicesList[j]
        ind1 = indices[0]
        c1 = c[ind1, statesKeep].toarray()[0] + unmerged[
                                         ind1] * unmerged[statesKeep] * 1.0 / c.shape[0]
        # BACE BF inverted so can use sparse matrices
        d[j, :indices[1].shape[0]] = 1. / multiDistSparseHelper(
            indices[1], c1, w[ind1], c, w, statesKeep, unmerged)

    return d


def multiDistSparseHelper(indices, c1, w1, c, w, statesKeep, unmerged):
    d = np.zeros(indices.shape[0], dtype=np.float32)
    p1 = c1 / w1
    for i in xrange(indices.shape[0]):
        ind2 = indices[i]
        c2 = c[ind2, statesKeep].toarray()[0] + unmerged[
                                         ind2] * unmerged[statesKeep] * 1.0 / c.shape[0]
        p2 = c2 / w[ind2]
        cp = c1 + c2
        cp /= (w1 + w[ind2])
        d[i] = c1.dot(np.log(p1 / cp)) + c2.dot(np.log(p2 / cp))
    return d

def calcSurprisalSparse(indices, cMats, w_multi, statesKeep, uncombined):

    d = 1./ surprisalSparseHelper(indices, cMats, w_multi, statesKeep, uncombined)

    return d


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
    preCombine = preCombineSparse
    multiDist = multiDistSparse
    calcSurprisal = calcSurprisalSparse
    run(cMats=cMats, nMacro=2, nProc=4, preCombine=preCombine, multiDist=multiDist, calcSurprisal=calcSurprisal)

if __name__ == "__main__":
    test_precombine()




