import os,sys,commands
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
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

cutoff = ct.Cutoff(3.0,8.0,0.5)
cutoffs = cutoff.ext(4.5,6,0.1)
cutoffs = cutoff.cutoffs

print cutoffs

metrics = 'dihedral'
if metrics.lower() == 'dihedral':
    Path = "/Users/tud51931/projects/MSM/msm/ff03-dihedralhybrid/"
    metrics = 'Dihedral'
    plotxlabel = '%s Cutoffs(Degrees)'%metrics
elif metrics.lower() == 'rmsd':
    Path = "/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/"
    metrics = 'RMSD'
    plotxlabel = '%s Cutoffs(Angstrom)'%metrics

if 1:
    
    numberofstates = []    
    
    for cutoff in cutoffs:
        metricdir = '%sCluster%0.1f'%(metrics,cutoff)
        targetDataDir = os.path.join(Path,metricdir,'lagtime50') 
        if not os.path.exists(targetDataDir):
            print "Can't find target directory:%s"%targetDataDir
            sys.exit()
        if os.path.exists(targetDataDir):
            print "Now I am in %s"%targetDataDir
            os.chdir(targetDataDir)
            if os.path.exists('Populations.dat'):
               output = commands.getoutput('wc Populations.dat')
               outputlines = output.split()
               numberofstates.append(int(outputlines[0]))
    print 'Cutoffs:',cutoffs
    print 'Number of States:',numberofstates
    
    if 1 :
        plt.figure()
        plt.xlabel(plotxlabel)
        plt.ylabel('Number of States')
        plt.plot(cutoffs,numberofstates)
        plt.title('Number of States-cutoffs')
        #plt.show()
        plt.savefig('numstates-%scutoffs.pdf'%metrics)
