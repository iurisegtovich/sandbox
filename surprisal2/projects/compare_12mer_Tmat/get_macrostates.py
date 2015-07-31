#! /usr/bin/env python

import os, sys
import numpy as np

from scipy.io import mmread, mmwrite

usage = """Usage: getMacrostates.py microstates.dat outdir
    
   OUTPUT
   outdir/macrostates.dat    - list of indices and conact states
   outdir/U_states.dat       - index of the unfolded contact state 
   outdir/F_states.dat       - index of the folded contact state """

sys.path.append('../../scripts')
from HelperTools  import *

################
# Main program

if len(sys.argv) < 3:
    print usage
    sys.exit(1)

microstatesFn = sys.argv[1]
outdir        = sys.argv[2]
if not os.path.exists(outdir):
    os.mkdir(outdir)

m = MicrostateInfo(microstatesFn)

outfile = os.path.join(outdir,'macrostates.dat')
fout = open(outfile,'w')
for i in range(len(m.uniqueContactStates)):
    cstate = m.uniqueContactStates[i]
    outstr = '%d\t%r'%(i,cstate)
    print outstr
    fout.write(outstr+'\n')
print 'Wrote %s -- Done.'%outfile 

# Find U and F states as max and min contacts
numcontacts = np.array([len(cstate) for cstate in m.uniqueContactStates])
U_state = numcontacts.argmin()
F_state = numcontacts.argmax()
print 'U_state', U_state, 'F_state', F_state
np.savetxt(os.path.join(outdir,'U_states.dat'), np.array([U_state]), fmt='%d')
np.savetxt(os.path.join(outdir,'F_states.dat'), np.array([F_state]), fmt='%d')




