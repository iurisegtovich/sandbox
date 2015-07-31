#!/usr/bin/env python

import os, sys
import numpy as np

from msmbuilder import MSMLib
from scipy.io import mmread
from scipy import loadtxt, savetxt

usage = """Usage:  python get_populations.py tProbsFn.mtx outdir
    will write tProbsFn.Populations.dat"""

if len(sys.argv) < 3:
    print usage
    sys.exit(1)   

tProbsFn = sys.argv[1]
outdir = sys.argv[2]
TC = mmread(tProbsFn)
EigAns = MSMLib.GetEigenvectors(TC,5)
Populations = EigAns[1][:,0]

outfile = os.path.join(outdir, os.path.basename(tProbsFn).replace('.mtx','.Populations.dat'))
print 'outfile', outfile
savetxt(outfile, Populations)

