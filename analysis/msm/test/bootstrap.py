import sys,os
import numpy as np
from msmbuilder import Serializer
from random import choice

"""
Input: A set of data to be boot strapped.

output: 
"""

def sample_with_replacement(datalist,numtraj):    

    return [choice(zip(datalist,range(len(datalist)))) for _ in range(numtraj)]

def bootstrap(AssignmentsFile,numtraj,bootstrapnumber,PathtoSaveFiles):
    bootstraplist,TrajID = [],[]
    File = Serializer.LoadFromHDF(AssignmentsFile)
    datalist = File['Data']
    if isinstance(numtraj,str) and numtraj.lower()=='all':
        numtraj = len(datalist)
    elif int(numtraj) <= 0:
        print "Please input valid number from 1 to %d"%len(datalist)
        sys.exit()
    else:
        numtraj = int(numtraj)
    
    for i in range(bootstrapnumber):
        all_assignments = -1 * np.ones((numtraj,len(datalist[0])), dtype=np.int)
        k = 0
        trajid = []
        for j,m in sample_with_replacement(datalist,numtraj):        
            all_assignments[k][:] = j[:]
            trajid.append(m)
            k += 1
        TrajID.append(trajid)
        bootstraplist.append(all_assignments)
    SaveBootstrapFiles(bootstraplist,TrajID,PathtoSaveFiles,File,AssignmentsFile)
        
def SaveBootstrapFiles(bootstraplist,TrajID,PathtoSaveFiles,File,AssignmentsFile):
    
    if not os.path.exists(PathtoSaveFiles):
        os.makedirs(PathtoSaveFiles)
    i = 0
    for data,trajid in zip(bootstraplist,TrajID):
        File['Data'] = data
        File['TrajID'] = trajid
        NewFile = AssignmentsFile.split('/')[-1].split('.')
        NewFile.insert(-1,'bs')
        NewFile[0] += '.%d'%i
        NewFile = '.'.join(NewFile)
        NewFile = os.path.join(PathtoSaveFiles,NewFile)
        try:
            File.SaveToHDF(NewFile)
            print "Save To :%s "%NewFile
        except Exception:
            print "File %s Already Exists! Pass it."
            pass
        i += 1    
    
def main():
    if len(sys.argv) != 5:
        print "Usage: python bootstrap.py Path to Assignments.h5 number of trajectory bootstrapnumber path to save files"
        sys.exit()
    bootstrap(AssignmentsFile=sys.argv[1],numtraj=sys.argv[2],bootstrapnumber=int(sys.argv[3]),PathtoSaveFiles=sys.argv[4])
if __name__=="__main__":
    main()



    
   
