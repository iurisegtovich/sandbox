import os,sys
import numpy as np
sys.path.append('/Users/tud51931/scripts/gfzhou')
import cutoffs as ct
RUN = True # If Ture, run the cmd, False, print cmd
DEBUG = False
TESTING = RUN

def run_cmd(cmd,testing=False):

    """ If true,RUN, run the cmd, if False print the cmd"""

    if testing:
        os.system(cmd)
    else:
        print '>>',cmd
        
def LinkData(sourceDataDir,targetDataDir):
    nruns = 100
    if not os.path.exists('%s/Trajectories'%sourceDataDir):
        print "Can't find the source -> Trajectories"
        sys.exit()        
        
    for run in range(nruns):
        sourcelh5 = os.path.join('%s/Trajectories'%sourceDataDir, 'trj%d.lh5'%run)
        if os.path.exists(sourcelh5):
            if not os.path.exists('%s/Trajectories'%targetDataDir):
                os.mkdir('%s/Trajectories'%targetDataDir)        
            if os.path.exists('%s/Trajectories'%targetDataDir):
                targetlh5 = os.path.join('%s/Trajectories'%targetDataDir, 'trj%d.lh5'%run)
                os.system('ln -s %s %s'%(sourcelh5, targetlh5))  # make a soft symbolic link                
                
    if not (os.path.exists('%s/AtomIndices.dat'%sourceDataDir) and os.path.exists('%s/gen0.pdb'%sourceDataDir) and os.path.exists('%s/ProjectInfo.h5'%sourceDataDir)):
        print "Can't find the source files -> AtomIndices.dat or gen0.pdb or ProjectInfo.h5"
        sys.exit()
    else:
        run_cmd('cp %s/AtomIndices.dat %s/gen0.pdb %s/ProjectInfo.h5 %s'%(sourceDataDir,sourceDataDir,sourceDataDir,targetDataDir),TESTING)


cutoff = ct.Cutoff(4.6,6.0,0.1)
#cutoffs = cutoff.ext(3.0,4.5,0.1)
cutoffs = cutoff.cutoffs
print cutoff.cutoffs

if 1:
    for cutoff in cutoffs:
        path = '/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter'
        sourceDataDir = '/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/sourcedata'
        targetDataDir = os.path.join(path,'RMSDCluster%0.1f'%cutoff) 
        if not os.path.exists(targetDataDir):
            os.mkdir(targetDataDir)
        if os.path.exists(targetDataDir):
            LinkData(sourceDataDir,targetDataDir)
