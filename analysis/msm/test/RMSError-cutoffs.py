import os,sys
sys.path.append('/Users/tud51931/scripts/gfzhou/')
from cutoffs import Cutoff
from RMSError import CalculateRMSError
import numpy as np
from scipy import loadtxt
import matplotlib.pyplot as plt
import commands


def convert_populations(Populations,Mapping):
    p = np.zeros(len(Mapping))
    for i in range(len(Mapping)):
        if int(Mapping[i]) != -1:
            try:
                p[i] = Populations[int(Mapping[i])]
            except IndexError:
                pass
    return p


tau = 50
cutoff = Cutoff(2.0,6.0,0.5)
#cutoffs = cutoff.ext(2.0,3.0,0.5)
cutoffs = cutoff.cutoffs
#cutoffs  =[2.0,3.0,3.5,4.0]
print 'cutoffs:',cutoffs

metrics = 'rmsd'
if metrics.lower() == 'dihedral':
    Path = "/Users/tud51931/projects/MSM/msm/ff03-dihedralhybrid/"
    metrics = 'Dihedral'
    plotxlabel = '%s Cutoffs(Degree)'%metrics
elif metrics.lower() == 'rmsd':
    Path = "/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/"
    metrics = 'RMSD'
    plotxlabel = '%s Cutoffs(Angstrom)'%metrics

plottitle = 'RMSError of State Populations - %s Cluster'%metrics
plotylabel = 'RMSError of State Populations'

E_ktau_cutoffs,E_ktau_rel_cutoffs = [],[]
RMSError_cutoffs,RMSError_mean_cutoffs,RMSError_std_cutoffs = [],[],[]
bs = False
if bs == False:
    for cutoff in cutoffs:
    
        metricdir = '%sCluster%0.1f'%(metrics,cutoff)
        print "Calculating %s now."%metricdir
        RMSError=[]
        filepath = os.path.join(Path,metricdir)
        p = os.path.join(filepath,'lagtime%d'%tau,'Populations.dat')
        m = os.path.join(filepath,'lagtime%d'%tau,'Mapping.dat')
        p_raw = os.path.join(filepath,'Data','Populations_raw.dat')
        try:
            Populations = loadtxt(p)
            Mapping = loadtxt(m)
        except IOError:
            print "Please build MSM first."
            raise
        try:
            Populations_raw = loadtxt(p_raw)
        except IOError:
            a = os.path.join(filepath,'Data','Assignments.h5')
            print "python ~/scripts/gfzhou/CalculateRawStatePopulations.py -a %s -o %s"%(a,p_raw)
            os.system("python ~/scripts/gfzhou/CalculateRawStatePopulations.py -a %s -o %s"%(a,p_raw))
            Populations_raw = loadtxt(p_raw)
        Populations = convert_populations(Populations,Mapping)
       
            
        RMSError_cutoffs.append(CalculateRMSError(Populations,Populations_raw))
      
    print 'RMSError:',RMSError_cutoffs
    plt.figure()
    plt.title(plottitle)
    plt.ylabel(plotylabel)
    plt.xlabel(plotxlabel)
    plt.plot(cutoffs,RMSError_cutoffs) 
    plt.show()
elif bs == True:
    for cutoff in cutoffs:
        metricdir = '%sCluster%0.1f'%(metrics,cutoff)
        print "Calculating %s now."%metricdir
        RMSError=[]
        filepath = os.path.join(Path,metricdir)
        path_bs = os.path.join(filepath,'lagtime%d'%tau,'bootstrap','bs0')
        pfile = commands.getoutput('ls %s/Populations*'%path_bs).split('\n')
        p_raw = os.path.join(filepath,'Data','Populations_raw.dat')
        for p in pfile:
            m = p.replace('Populations','Mapping')
            try:
                Populations = loadtxt(p)
                Mapping = loadtxt(m)
            except IOError:
                print "Please build MSM first."
                raise
            try:
                Populations_raw = loadtxt(p_raw)
            except IOError:
                a = os.path.join(filepath,'Data','Assignments.h5')
                print "python ~/scripts/gfzhou/CalculateRawStatePopulations.py -a %s -o %s"%(a,p_raw)
                os.system("python ~/scripts/gfzhou/CalculateRawStatePopulations.py -a %s -o %s"%(a,p_raw))
                Populations_raw = loadtxt(p_raw)
            
            Populations = convert_populations(Populations,Mapping)
            Mapping = loadtxt(os.path.join(filepath,'lagtime%d'%tau,'Mapping.dat'))
            Populations = convert_populations(Populations,Mapping)
            RMSError.append(CalculateRMSError(Populations,Populations_raw))
        RMSError_mean_cutoffs.append(np.mean(RMSError))
        RMSError_std_cutoffs.append(np.std(RMSError))
      
    print 'RMSError_mean:',RMSError_mean_cutoffs
    print 'RMSError_std', RMSError_std_cutoffs
    plt.figure()
    plt.title(plottitle)
    plt.ylabel(plotylabel)
    plt.xlabel(plotxlabel)
    plt.errorbar(cutoffs,RMSError_mean_cutoffs,RMSError_std_cutoffs) 
    plt.savefig('RMSError-cutoffs_bs.png')
    plt.show()    
    
