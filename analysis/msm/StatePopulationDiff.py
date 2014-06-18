#-------------------------
#This script is to calculate the State Population Difference between state 
#population from MSM and state population from raw data.
#The result will be a PopulationDifference vs cutoffs plot. 
#To some extent it will show how kinetic accurate and 
#statistical reliable is for each model with different cutoffs.
#-------------------------
#CHANGE LOG
#Changed the script name from RMSE.py to StatePopulationDiff.py
#Guangfeng,March,5,2013
#-------------------------
#Guangfeng Zhou
#Voelz Lab
#Chemistry Department
#Temple University
#-------------------------
"""
TODO: Add another option argument to allow user defined cutoff list and directory list.
"""

import os,sys
import argparse
import commands
import numpy as np
from scipy import loadtxt
import matplotlib.pyplot as plt


def CalculateRMSE(inputdir,metric,lagtime):
    RMSError_cutoffs = []
    tau = lagtime
    Dirs =  commands.getoutput('ls %s | grep -i %s'%(inputdir,metric)).split('\n')
    cutoffs = getcutoffs(Dirs)
    for (cutoff,Dir) in zip(cutoffs,Dirs):
        print "Calculating %s now."%Dir
        RMSError=[]
        filepath = os.path.join(inputdir,Dir)
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
            print "Can't find Populations_raw.dat in Data/ directory.Please Run CalculateRawPopulations.py first to get Populations_raw.dat."
            print "python ~/scripts/gfzhou/CalculateRawStatePopulations.py -a %s -o %s"%(a,p_raw)
            os.system("python ~/scripts/gfzhou/CalculateRawStatePopulations.py -a %s -o %s"%(a,p_raw))
            Populations_raw = loadtxt(p_raw)
        Populations = convert_populations(Populations,Mapping)
        RMSError_cutoffs.append(CalculateRMSError(Populations,Populations_raw))        
    return cutoffs,RMSError_cutoffs

def CalculateRMSError(Populations,Populations_raw):
     
    p = loadtxt(Populations)
    p_raw = loadtxt(Populations_raw)
    if len(p) != len(p_raw):
        print "Populations.dat doesn't match Populations_raw.dat"
        sys.exit()
    else:
        rmserror = (sum([(p[i]-p_raw[i])**2 for i in range(len(p))])/len(p))**0.5 # Note: It is not divided by the total number of states
    return rmserror

def convert_populations(Populations,Mapping):
    p = np.zeros(len(Mapping))
    for i in range(len(Mapping)):
        if int(Mapping[i]) != -1:
            try:
                p[i] = Populations[int(Mapping[i])]
            except IndexError:
                pass
    return p

def getcutoffs(Dirs):
    """
    Get a list of cutoffs from directories. Only works when the directories have certain names ending with cutoffs such as RMSDCluster4.2, DihedralCluster3.0
    """
    cutoff = 0
    cutoffs = []
    for Dir in Dirs:
        for char in Dir:
            if char.isdigit():
                cutoff = float(Dir[Dir.find(char):])
                break
        cutoffs.append(cutoff)
    return cutoffs

def Function():
    parser = argparse.ArgumentParser()
    parser.add_argument('-m','--Metric',
                        help="The metric used to cluster data.",
                        type=str,choices=['rmsd','dihedral'],
                        required=True)
    parser.add_argument('-l','--Lagtime',
                        help="The Lagtime used to build MSM.",
                        type=int,default=50)
    parser.add_argument('-i','--Input_Dir',
                        help="Path to the parent directory containing subdirectories which contains Markov State Model directories.",required = True)
    args = parser.parse_args()
    
    cutoffs,RMSError_cutoffs=CalculateRMSE(args.Input_Dir,args.Metric,args.Lagtime)
    
    print 'RMSError:',RMSError_cutoffs
    
    plotylabel = 'RMSError of State Populations'  
    plotxlabel = '%s Cutoffs(Angstrom)'%args.Metric.upper()
    plottitle = 'RMSError of State Populations - %s Cluster'%args.Metric.upper()
    plt.figure()
    plt.title(plottitle)
    plt.ylabel(plotylabel)
    plt.xlabel(plotxlabel)
    plt.plot(cutoffs,RMSError_cutoffs) 
    #plt.show()    
    plt.savefig('RMSErrorofStatePopulation-FsERR.png')

if __name__ == "__main__":
    Function()
