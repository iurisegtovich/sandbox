import os
import commands
import numpy as np
import argparse
from scipy.io import mmread
from scipy import loadtxt,savetxt

parser = argparse.ArgumentParser()
parser.add_argument('-i','--Input_Dir',help="Path to directory containing all tCounts.mtx files. Note: ONLY .mtx files")
parser.add_argument('-o','--Output_Dir',help="Output directory containing two population files. One is state equilibrium populations over all tCounts file, the other one is the average state populations calculated from the first file. Default: Input_Dir/Population/")

args = parser.parse_args()

fns = commands.getoutput('ls %s'%args.Input_Dir).split('\n')
for tcount in fns:
    fn = os.path.join(args.Input_Dir,tcount)
    c = mmread(fn)
    p = c.sum(1)/c.sum()
    try:
        P = np.vstack((P,p.reshape(1,-1)))
    except NameError,ValueError:
        P = p.reshape(1,-1)

if args.Output_Dir == None:
    args.Output_Dir = os.path.join(args.Input_Dir,'Population')
if not os.path.exists(args.Output_Dir):
    os.makedirs(args.Output_Dir)
pops = os.path.join(args.Output_Dir,'Populations.dat')
meanpops = os.path.join(args.Output_Dir,'MeanPopulations.dat')
savetxt(pops,P)
print "Wrote %s"%pops
savetxt(meanpops,P.mean(0))
print "Wrote %s"%meanpops
