# useSurprisal.py

import os,sys
sys.path.append('/Users/tud51931/scripts/gfzhou/')
import Surprisal
import commands
import argparse
import numpy as np
from scipy import savetxt,loadtxt
import matplotlib.pyplot as plt
from scipy.io import mmread

parser = argparse.ArgumentParser()
parser.add_argument('-oc','--OCount',help="Original UnMapped tCount.mtx file")
parser.add_argument('-bc','--BCount',help="Path to bootstrapped tCount.mtx file",action = 'append')
parser.add_argument('-c','--cutoff',help="help to identify which data are used",type=float,required=True)
parser.add_argument('-m','--metric',help="help to identify which data are used",type=str,required=True)
parser.add_argument('-n','--numtrajs',help="the number of trajectories used",type = int,required=True)
parser.add_argument('-s','--savesurprisal',action = 'store_true',default=True)
parser.add_argument('-p','--populations',help="Populations.dat file",required = True)
parser.add_argument('-d','--dvalue',action = 'store_true',default=True)

args = parser.parse_args()

sparse1 = mmread(args.OCount)

bcounts = []
for path in args.BCount:
    if path != '':
        mtxfiles = os.path.join(path,'*.mtx')
        fns = commands.getoutput('ls %s'%mtxfiles).split('\n')
        bcounts += fns
numoftcounts = len(bcounts)
surprisals_bootstraps = np.zeros((numoftcounts,sparse1.shape[0]))

for i in range(numoftcounts):
    print "Calculate %d of %d tCounts"%(i,numoftcounts)
    sparse2 = mmread(bcounts[i])
    obj = Surprisal.SurprisalCalculator(sparse1, sparse2)
    surprisals = obj.calculate_all_surprisal(weighted=True)
    surprisals_bootstraps[i,:] = surprisals[:]
surprisals_bootstraps_mean = surprisals_bootstraps.mean(0)
print len(surprisals_bootstraps)
surprisals_bootstraps_std = surprisals_bootstraps.std(0)
if args.dvalue:
    pop = loadtxt(args.populations)
    dvalue = (pop*surprisals_bootstraps).sum(1)
    fname = 'dvalues%s%0.1f.%dtrajs.dat'%(args.metric,args.cutoff,args.numtrajs)
    savetxt(fname,dvalue)
    print "Save to %s"%fname
print "dvalue",dvalue
"""
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
"""
