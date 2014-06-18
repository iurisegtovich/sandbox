#getUnmappedCount.py
#use GetCount.py written by Brandon

import os,sys
import commands
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-a','--Assignments',help='''Path to the folder of Assignment files. Attention: this should be a path''')
parser.add_argument('-l','--lagtime',help='''Lagtime,default = 50 (steps)''',default=50,type=int)
parser.add_argument('-o','--output',help='''Path to the output directory,default = .yourPathToAssignemts/UnMappedCounts''')

args = parser.parse_args()
args.output = os.path.join(args.Assignments,'UnMappedCounts')
if not os.path.exists(args.output):
    os.makedirs(args.output)
assignmentsfiles = os.path.join(args.Assignments,'*.h5')
assignments = commands.getoutput('ls %s'%assignmentsfiles).split('\n')

for a in assignments:
    afn = a.split('/')[-1]
    cfn = afn.replace('Assignments','tCounts')
    cfn = cfn.replace('h5','mtx')
    output = os.path.join(args.output,cfn)
    os.system('python ~/voelzlab/sandbox/gfzhou/test/GetCount.py %s %d %s'%(a,args.lagtime,output))
    print "Save to:%s"%output
    
