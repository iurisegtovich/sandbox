#-------------------------------------
#CalculateRawStatePopulations.py
#Written by Guangfeng Zhou
#Dr.Voelz Lab
#Chemistry Department
#Temple University
#November,1,2012
#--------------------------------------

import os,sys
import argparse
from msmbuilder import Serializer
from scipy import savetxt
import numpy as np

def calculate_statepopulation_rawdata(AssignmentsFixed):
    
    a = Serializer.LoadFromHDF(AssignmentsFixed)
    statenumber = max([max(a['Data'][i]) for i in range(len(a['Data']))])+1
    p = np.zeros(statenumber)
    for state in range(statenumber):
        for traj in range(len(a['Data'])):
            p[state] += a['Data'][traj].tolist().count(state)
    p = p/p.sum()
    return p

parser = argparse.ArgumentParser()
parser.add_argument('-a','--Assignment',help="Path to Assigenments File,default = Data/Assignments.h5",default = "Data/Assignments.h5")
parser.add_argument('-o','--Output',help="Output file, default=Populations_raw.dat",default="Populations_raw.dat")
args = parser.parse_args()
StatePopulation_rawdata = calculate_statepopulation_rawdata(args.Assignment)
savetxt(args.Output,StatePopulation_rawdata)
print "Save to %s"%args.Output
print "Done."
