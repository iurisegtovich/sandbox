import os
from scipy.io import mmread
import argparse
import commands

parser = argparse.ArgumentParser()
parser.add_argument('-i','--Input_Dir',help="Directory contains .mtx files",required=True)
args = parser.parse_args()

path = os.path.join(args.Input_Dir,'*.mtx')
fns = commands.getoutput('ls %s'%path).split('\n')
for fn in fns:
    M = mmread(fn)
    print fn,M.shape