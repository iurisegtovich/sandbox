import os, sys, glob

import numpy as np
import scipy.io
from scipy import loadtxt, savetxt

sys.path.append('../src')
from Surprisal import *

from matplotlib import pyplot as plt

def most_lumpable_pair(c):
    """Computes the surprisal for all pairs, and returns the pair with the smallest value

    INPUT
    c		an array C_ij of transition counts i -> j

    RETURNS
    (i,j)  	the pair with the smallest surprisal value
    S           the (unnormalized) surprisal value,
    var_S       ...and it's variance
    """
    # Find the most lumpable pair
    best_S = 1.e99
    best_pair = (None, None)
    for i in range(c.shape[0]): 
      for j in range(i+1, c.shape[0]):  
        S = surprisal(c[i,:], c[j,:], normalize=False)
        #print 'S_%d,%d ='%(i,j), S,
        if S < best_S:
            best_pair = (i,j)
            best_S = S
            #print '<-- best so far'
        #else:
        #    print

    # estimate the variance of S
    c1 = c[best_pair[0],:].astype(int)
    c2 = c[best_pair[1],:].astype(int)
    # print 'c1', c1, 'c2', c2
    var_S = surprisal_var(c[best_pair[0],:], c[best_pair[1],:], normalize=False)

    return best_pair, best_S, var_S
    
    
def lump(c, i, j):
    """For a count matrix c, lump states i and j together.

    RETURNS
    c_lumped		the lumped count matrix.

    WARNING
    This only works for dense matrices for now (June 2014)
    """

    n = c.shape[0]

    # Order indices so j is the last row and column
    Ind = range(n)
    Ind.pop(j)
    Ind.append(j)
    d = c[Ind,:]  # reorder the rows
    d = d[:,Ind]  # reorder the columns

    # add the last (was jth) row and column to the ith
    d[i,:] += d[-1,:]
    d[:,i] += d[:,-1]

    # return the lumped counts as a n-1 x n-1 submatrix
    return d[0:n-1,0:n-1]

def neglogP(C):
    """Returns f = -ln P(C|M) and the estimate of the variance var(L)."""

    # get the number of states
    n = C.shape[0] 

    # -ln P(C|M) = \sum_i -ln P(C_i|M) 
    f = 0.0
    f_var = 0.0
    for i in range(n):
        C_i = C[i,:].sum()
        f += C_i*H(C[i,:]) #+ n*np.log(n) - n
        f_var += C_i*H_var(C[i,:])

    return f, f_var

nmacro = 20  # 200
tCountFn = 'wt-bb-cb-2000-2ev-t-%dmacro-tCounts.mtx'%nmacro
#tCountFn = 'wt-bb-cb-2ev-t-2000micro-tCounts.mtx'    # 2000 states!
c_wt = np.array( scipy.io.mmread(tCountFn).todense() ) 

tCountFn = 'tz4-bb-cb-2000-2ev-t-%dmacro-tCounts.mtx'%nmacro
#tCountFn = 'tz4-bb-cb-2ev-t-2000micro-tCounts.mtx'   # 2000 states!
c_mut = np.array( scipy.io.mmread(tCountFn).todense() )

# Keep lumping according to BACE
c = c_wt + 1.0   # VAV: This is to compare with BACE, which uses pseudocounts
L = [[i] for i in range(c.shape[0])]   # keep track of lumping [[0], [1,4], [2], [3], [5,6,7]]
neglogBF = [0.0] 
variances = [0.0]
neglogPs = []
neglogPs_var = []
while c.shape[0] > 1:
    b, b_var = neglogP(c)
    print '*** -ln P(C|M) =', b, '+/-', b_var
    neglogPs.insert(0,b)
    neglogPs_var.insert(0,b_var)
    (i, j), best_S, var_S = most_lumpable_pair(c)
    neglogBF = [best_S] + neglogBF
    variances = [var_S] + variances
    print c.shape[0], (i, j), 'S =', best_S, '+/-', var_S  
    c = lump(c, i, j)
    L[i] += L.pop(j)
    #print 'Lumping', L
print 'neglogBF',neglogBF,
plt.figure()

plt.subplot(2,2,1)
#plt.plot(range(1,len(neglogBF)+1), neglogBF, '-')
plt.errorbar(range(1,len(neglogBF)+1), neglogBF, yerr=variances)
plt.hold(True)
bace_data = loadtxt('Output_BACE/bayesFactors.dat')
plt.plot(bace_data[:,0]-1, bace_data[:,1], 'r-')

plt.subplot(2,2,2)
plt.errorbar(range(1,len(neglogPs)+1), neglogPs, yerr=neglogPs_var)

plt.subplot(2,2,4)
plt.plot(range(1,len(neglogPs_var)+1), neglogPs_var)
plt.show()

print 'For lumping WT with MUT:'
for i in range(20): 
    print 'S_%d ='%i, surprisal(c_wt[i,:], c_mut[i,:], normalize=False)



