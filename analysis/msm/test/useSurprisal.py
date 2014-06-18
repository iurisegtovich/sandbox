# useSurprisal.py

import os,sys
sys.path.append('/Users/tud51931/scripts/gfzhou/')
import Surprisal
import commands
import argparse
import numpy as np
from scipy import savetxt
import matplotlib.pyplot as plt
from scipy.io import mmread

parser = argparse.ArgumentParser()
parser.add_argument('-oc','--OCount',help="Original UnMapped tCount.mtx file")
parser.add_argument('-bc','--BCount',help="Path to bootstrapped tCount.mtx file",action = 'append')
parser.add_argument('-c','--cutoff',help="help to identify which data are used",type=float)
parser.add_argument('-m','--metric',help="help to identify which data are used",type=str)
parser.add_argument('-n','--numtrajs',help="the number of trajectories used",type = int)
parser.add_argument('-s','--savesurprisal',action = 'store_true')

args = parser.parse_args()

sparse1 = mmread(args.OCount)

bcounts = []
for path in args.BCount:
    if path != '':
        fns = os.path.join(path,'*.mtx')
        fns = commands.getoutput('ls %s'%fns).split('\n')
        bcounts += fns
numoftcounts = len(bcounts)
surprisals_bootstraps = np.zeros((numoftcounts,sparse1.shape[0]))

for i in range(numoftcounts):
    print "Calculate %d of %d tCounts"%(i,numoftcounts)
    sparse2 = mmread(bcounts[i])
    obj = Surprisal.SurprisalCalculator(sparse1, sparse2)
    surprisals = obj.calculate_all_surprisal(weighted=True)
    surprisals_bootstraps[i,:] = surprisals[:]
surprisals_bootstraps = surprisals_bootstraps.transpose()
surprisals_bootstraps_mean = surprisals_bootstraps.mean(1)
surprisals_bootstraps_std = surprisals_bootstraps.std(1)
print "mean",surprisals_bootstraps_mean
print "std",surprisals_bootstraps_std
plt.figure()
plt.errorbar(range(len(surprisals_bootstraps_mean)),surprisals_bootstraps_mean,surprisals_bootstraps_std,fmt='r.')
plt.xlabel('State')
plt.ylabel('Surprisal Value')
plt.title('Surprisal Value - State')
figname = 'SurprisalValue_State_%dtcounts%s%0.1f.%dtrajs.weighted.png'%(numoftcounts,args.metric,args.cutoff,args.numtrajs)
plt.savefig(figname)
print "Save to %s"%figname
if args.savesurprisal:
    fn = figname.replace('png','dat')
    savetxt(fn,surprisals_bootstraps_mean)
    print "Save to %s"%fn
#plt.show()
