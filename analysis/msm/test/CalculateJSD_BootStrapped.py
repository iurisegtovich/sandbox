# CalculateJSD_BootStrapped.py

import os,sys
sys.path.append('/Users/tud51931/voelzlab/surprisal/')
from Surprisal import SurprisalCalculator, get_populations
import commands
import argparse
import numpy as np
from scipy import savetxt,loadtxt
import matplotlib.pyplot as plt
from scipy.io import mmread

parser = argparse.ArgumentParser()
parser.add_argument('-oc','--OCount',help="Original UnMapped tCount.mtx file",required=True)
parser.add_argument('-bc','--BCount',help="Path to bootstrapped tCount.mtx file",action = 'append',required=True)
parser.add_argument('-m','--Popmethod',choices=['mle','counts','retain','retain_first_mle'],help="Method to calculate populations.Default:'mle'",default = 'mle')
parser.add_argument('-f','--Fix', action='store_true', help='Whether to fix sparse matrix.Default=True',default=True)
parser.add_argument('-o','--Output',help="Output file.Default:JSDs.dat",default='JSDs.dat')

args = parser.parse_args()

sparse1 = mmread(args.OCount)

bcounts = []
for path in args.BCount:
    if path != '':
        mtxfiles = os.path.join(path,'*.mtx')
        fns = commands.getoutput('ls %s'%mtxfiles).split('\n')
        bcounts += fns
numoftcounts = len(bcounts)
JSDs_bootstraps = np.zeros((numoftcounts,sparse1.shape[0]))

if args.Popmethod == 'retain_first_mle':
    Populations1 = get_populations(count_mtx = sparse1, method = 'mle')[0]
    #for testing
    #Populations1 = get_populations(count_mtx = sparse1, method = 'counts')
else:
    Populations1 = None

for i in range(numoftcounts):
    print "Calculate %d of %d tCounts"%(i,numoftcounts)
    sparse2 = mmread(bcounts[i])
    print sparse1.shape,sparse2.shape
    if args.Fix == True:
        fix = mmread(args.OCount)
        fix = fix/fix.sum()
        fix = fix.asformat('csr')
        for j in range(fix.shape[0]):
            for k in range(fix.shape[0]):
                if fix[j,k] != 0:
                    fix[j,k] = np.ceil(fix[j,k])
        fix = fix.asformat('coo')
        sparse2 = sparse2 + fix
    obj = SurprisalCalculator(sparse1, sparse2, pop_method = args.Popmethod, init_guess1 = Populations1)
    JSDs = obj.calculate_all_jsd()
    JSDs_bootstraps[i,:] = JSDs[:]
print len(JSDs_bootstraps)
savetxt(args.Output,JSDs_bootstraps)
print "Wrote:%s"%args.Output
