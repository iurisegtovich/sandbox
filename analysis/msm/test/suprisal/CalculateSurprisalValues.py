# CalculateSurprisalValues.py
# For bootstraped data.

import os,sys
sys.path.append('/Users/tud51931/scripts/gfzhou/')
import Surprisal
import commands
import argparse
import numpy as np
from scipy import savetxt
from scipy.io import mmread

parser = argparse.ArgumentParser()
parser.add_argument('-oc','--OCount',help="Original UnMapped tCount.mtx file")
parser.add_argument('-bc','--BCount',help="Path to bootstrapped tCount.mtx file",action = 'append')
parser.add_argument('-w','--Weighted',action = "store_true")
parser.add_argument('-o','--Output',help="Output file containing all the surprisal values.Default: surprisalvalues.dat",default="surprisalvalues.dat")

args = parser.parse_args()

bcounts = []
sparse1 = mmread(args.OCount)
for path in args.BCount:
    if path != '':
        fns = os.path.join(path,'*.mtx')
        fns = commands.getoutput('ls %s'%fns).split('\n')
        bcounts += fns
numoftcounts = len(bcounts)
surprisals_bootstraps = np.zeros((numoftcounts,sparse1.shape[0]))
if args.Weighted == True:
    print "True, continue"
else:
    print "False, change"
    sys.exit()

for i in range(numoftcounts):
    print "Calculate %d of %d tCounts"%(i,numoftcounts)
    sparse2 = mmread(bcounts[i])
    obj = Surprisal.SurprisalCalculator(sparse1, sparse2)
    surprisals = obj.calculate_all_surprisal(weighted=args.Weighted)
    surprisals_bootstraps[i,:] = surprisals[:]

savetxt(args.Output,surprisals_bootstraps)
print "Wrote: %s"%args.Output