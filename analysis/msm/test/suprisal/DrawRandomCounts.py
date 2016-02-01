#DrawRandomCounts.py
#Only support one state now.
#-------------------
#Future endeavor:
#Pseudo sampling from multiple starting states with various distributions.(State Population based;AA;)
#--------------------
#Guangfeng Zhou
#Voelz Lab
#Chemistry Department, Temple University
#5/18/2013

import os,sys
import argparse
import random
import copy
import numpy as np
from msmbuilder import Serializer,Trajectory


def get_startpoint(stateid,assignment,frames):
    #print assignment['Data']
    rowindx,colindx = np.where(assignment['Data']==stateid)
    rowindx = rowindx.tolist()
    colindx = colindx.tolist()
    #random choose a valid starting point
    while True:
        if len(rowindx)==0:
            return None
        i = random.randint(0,len(rowindx)-1)
        if args.Verbose:
            print len(rowindx),i
        startpoint=(rowindx[i],colindx[i])
        try:
            if assignment['Data'][startpoint[0],startpoint[1]+frames] >=0:
                return startpoint
        except IndexError:
            if args.Verbose:
                print "Index Out Of Bounds"
            rowindx[i:i+1]=[]
            colindx[i:i+1]=[]
            pass
        
def pseudotrajs(numberoftrajs,stateid,assignment,frames):
    a = copy.copy(assignment)
    a.clear()
    for traj in range(args.NumberofTrajs):
        startpoint = get_startpoint(stateid,assignment,frames)
        if startpoint != None:
            try:
                a['Data'] = np.vstack((a['Data'],assignment['Data'][startpoint[0]][startpoint[1]:startpoint[1]+frames]))
                if args.Verbose:
                    print "vstack"
            except KeyError:
                a['Data'] = assignment['Data'][startpoint[0]][startpoint[1]:startpoint[1]+frames]
        else:
            return None
    return a

def pseudosampling(states,assignment,numberoftrajs,frames,output):
    try:
        fn = Serializer.LoadFromHDF(assignment)
    except IOError:
        print "Can't find Assignment file"  
        sys.exit()
    for stateid in states:
        a = pseudotrajs(numberoftrajs,stateid,fn,frames)
        a.SaveToHDF(output)
        print "Wrote:%s"%output



    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-s','--States',help="A list of surprising states",required=True,nargs='+',type=int)
    parser.add_argument('-a','--Assignments',help="Path to Assignments.h5 file",required=True)
    parser.add_argument('-n','--NumberofTrajs',help="The number of new sampling trajectories for each state. Default=1",default=1,type=int)
    parser.add_argument('-f','--Frames',help="The number of frames in each of the new trajecotries.Default = 50 (5ns)",default=50,type=int)
    parser.add_argument('-o','--Output',help="Pseudo Trajectory Assignment file.",required=True,type=str)
    parser.add_argument('-v','--Verbose',help="increase output verbosity",action="store_true")
    args = parser.parse_args()    
    
    pseudosampling(args.States,args.Assignments,args.NumberofTrajs,args.Frames,args.Output)
    