
import os,sys
sys.path.append('/Users/tud51931/scripts/gfzhou')
from cutoffs import Cutoff
import commands
cutoff = Cutoff(3.0,6.0,0.5)
#cutoffs = cutoff.ext(2.0,3.0,0.5)
cutoffs = cutoff.cutoffs
#cutoffs = [1.5]
print 'cutoffs:',cutoffs

for cutoff in cutoffs:
    for i in range(0,1):
        Path = "/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster%0.1f/lagtime50/bootstrap/bs%d"%(cutoff,i)
        os.chdir(Path)
        Assignmentsfile = commands.getoutput('ls *.h5').split('\n')
        for fn in Assignmentsfile:
            assignments = os.path.join(Path,fn)
            outputdir = os.path.join(Path,fn.replace('.Fixed.bs.h5',''))
            
            print outputdir
            #os.system("python ~/scripts/gfzhou/buildmsm_modified.py -l 50 -a %s --assignfilename %s"%(assignments,assignments))
            cmd = "BuildMSM.py -l 50 -a %s -o %s"%(assignments,outputdir)
            print cmd
            os.system(cmd)
