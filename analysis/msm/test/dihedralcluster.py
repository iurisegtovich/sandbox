import os,sys
import numpy as np
sys.path.append('Users/tud51931/scripts/gfzhou/')
from cutoffs import Cutoff

RUN = True # If Ture, run the cmd, False, print cmd
DEBUG = False
TESTING = RUN

def run_cmd(cmd,testing=False):

    """ If true,RUN, run the cmd, if False print the cmd"""

    if testing:
        os.system(cmd)
    else:
        print '>>',cmd



cutoff = Cutoff(4.5,6,0.1)
cutoffs = cutoff.cutoffs

print "cutoffs",cutoffs

BootStrap = False

for cutoff in cutoffs:
    path = '/Users/tud51931/projects/MSM/msm/ff03-dihedralhybrid'
    targetDataDir = os.path.join(path,'DihedralCluster%0.1f'%cutoff) 
    if not os.path.exists(targetDataDir):
        print "Can't find target directory:%s"%targetDataDir
        sys.exit()
    if os.path.exists(targetDataDir):
        print "Now I am in %s"%targetDataDir
        os.chdir(targetDataDir)
        #if not os.path.exists('../Plottimescale/cutoff%0.1f'%cutoff):
        #    os.mkdir('../Plottimescale/cutoff%0.1f'%cutoff)
        #else:
        #    run_cmd('cp ImpliedTimescales.dat AtomIndices.dat ProjectInfo.h5 ../Plottimescale/cutoff%s'%cutoff,TESTING)
        #run_cmd('Cluster.py -S 10 dihedral hybrid -d %3.3f'%(float(cutoff)))
        #run_cmd('Cluster.py -S 10 dihedral hybrid -d %3.3f'%(float(cutoff)),TESTING)
        #run_cmd('Assign.py dihedral')
        #run_cmd('Assign.py dihedral',TESTING)
        if BootStrap:
            for i in range(1,101):
                run_cmd("CalculateImpliedTimescales.py -a ./Data/assignments/Assignments%d.h5 -l 1,200 -i 5 -p 8 -o impliedtimescales/ImpliedTimescales%d.dat"%(i,i))
                run_cmd("CalculateImpliedTimescales.py -a ./Data/assignments/Assignments%d.h5 -l 1,200 -i 5 -p 8 -o impliedtimescales/ImpliedTimescales%d.dat"%(i,i),TESTING)
        elif 1:
            run_cmd("CalculateImpliedTimescales.py -l 1,200 -i 5 -p 8")
            run_cmd("CalculateImpliedTimescales.py -l 1,200 -i 5 -p 8",TESTING)
        #run_cmd('PlotImpliedTimescales.py',TESTING)
        #lagtimes = [150,200,250,300,350,400,500]            
        lagtimes = [50,100,150,200,250,300,350,400,450,500]
        for lagtime in lagtimes:
            if BootStrap:            
                for i in range(1,101):
                    #if os.path.exists('lagtime%d'%lagtime):
                    #    run_cmd('rm -rf lagtime%d'%lagtime,TESTING)
                    run_cmd('BuildMSM.py -a ./Data/assignments/Assignments%d.h5 -l %d -o lagtime%d/bootstrap/%d'%(i,lagtime,lagtime,i))
                    run_cmd('BuildMSM.py -a ./Data/assignments/Assignments%d.h5 -l %d -o lagtime%d/bootstrap/%d'%(i,lagtime,lagtime,i),TESTING)
            elif 0:
                run_cmd('BuildMSM.py -l %d -o lagtime%d'%(lagtime,lagtime))
                run_cmd('BuildMSM.py -l %d -o lagtime%d'%(lagtime,lagtime),TESTING)

