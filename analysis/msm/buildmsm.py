#----------------------
#Guangfeng Zhou
#Voelz Lab
#Chemistry Department
#Temple University
#---------------------


import os,sys
import argparse
import numpy as np

RUN = True # If Ture, run the cmd, False, print cmd
DEBUG = False
TESTING = RUN

def run_cmd(cmd,testing=False):

    """ If true,RUN, run the cmd, if False print the cmd"""

    if testing:
        os.system(cmd)
    else:
        print '>>',cmd

def buildmsm(cutoffs,targetDataDir,metric,lagtimes,symmetrize):
    for cutoff in cutoffs:
        if not os.path.exists(targetDataDir):
            print "Can't find target directory:%s"%targetDataDir
            sys.exit()
        if os.path.exists(targetDataDir):
            print "Now I am in %s"%targetDataDir
            os.chdir(targetDataDir)
            run_cmd('Cluster.py -S 10 %s -d %3.3f'%(metric,cutoff/10))
            run_cmd('Cluster.py -S 10 %s -d %3.3f'%(metric,cutoff/10),TESTING)
            run_cmd('Assign.py %s'%metric.split()[0])
            run_cmd('Assign.py %s'%metric.split()[0],TESTING)
            run_cmd('CalculateImpliedTimescales.py -l 1,200 -i 5 -p 8')
            run_cmd('CalculateImpliedTimescales.py -l 1,200 -i 5 -p 8', TESTING)
            for lagtime in lagtimes:
                run_cmd('BuildMSM.py -l %d -o lagtime%d -s %s'%(lagtime,lagtime,symmetrize))
                run_cmd('BuildMSM.py -l %d -o lagtime%d -s %s'%(lagtime,lagtime,symmetrize),TESTING)
    
def formydata(cutoffs,inputdir,metric,lagtimes,symmetrize):
    for cutoff in cutoffs:
        path = inputdir
        if metric.split()[0].lower() == 'rmsd':
            targetDataDir = os.path.join(path,'RMSDCluster%0.1f'%cutoff) 
        elif metric.split()[0].lower() == 'dihedral':
            targetDataDir = os.path.join(path,'DihedralCluster%0.1f'%cutoff)
        buildmsm(cutoff,targetDataDir,metric,lagtimes,symmetrize)

    
parser = argparse.ArgumentParser()
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-t','--TargetDir',
                   help="Target directory contains Trajectory folder, AtomIndices.dat, ProjectInfo.h5, pdb file")
group.add_argument('-i','--Input_Dir',
                    help="Path to the parent directory containing subdirectories which contains Trajectory folder, AtomIndices.dat, ProjectInfo.h5, pdb file. To use this option, certain folder structure is required.")
parser.add_argument('-c','--Cutoffs',nargs='*',type = float,
                    help="Specify the cutoff you want to use to cluster your data. Unit here is Angstrom. You can specify multiple cutoffs.",required = True)
parser.add_argument('-l','--Lagtimes',nargs='*',required = True,type=int,
                    help="Specify the lagtime you want to use to build Markov State Model. Unit here is step, which means if the time interval in your .xtc or .trr files is 100ps, then that lagtime = 50 means you use 50*100ps = 5ns to build MSM.Multiple lagimes can be specified.")
parser.add_argument('-s','--symmetrize',choices=['MLE','Transpose','None'],default = 'MLE',
                    help="Method by which to estimate a symmetric counts matrix.Symmetrization ensures reversibility, but may skew dynamics. We recommend maximum likelihood estimation(MLE) when tractable, else try Transpose. It isstrongly recommended you read the documentation surrounding this choice. Default: MLE")
parser.add_argument('-m','--Metric',default ='rmsd hybrid',type = str,
                    help="the method to cluster your data, please Run Cluster.py -h for further information. Default='rmsd hybrid'",required = True)
args = parser.parse_args()

if __name__=="__main__":
    if args.Input_Dir:
        formydata(args.Cutoffs,args.Input_Dir,args.Metric,args.Lagtimes,args.symmetrize)
    elif args.TargetDir:
        buildmsm(args.Cutoffs,args.TargetDir,args.Metric,args.Lagtimes,args.symmetrize)