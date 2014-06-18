#-------------------------------------------------
#Function:Reference_Rg was copied from msmbuilder 
#The rest was written by Guangfeng Zhou
#
#------CHANGE LOG-------
#Oct,24,2012, Fixed a bug and added the header, Guangfeng Zhou
#-----------------------
#
#Guangfeng Zhou
#Dr.Voelz Lab
#Chemistry Department
#Temple University
#October,24,2012


import sys,os
import numpy as np
from scipy import savetxt,loadtxt
from msmbuilder import Trajectory,Serializer
import argparse

def Reference_Rg(trajfile):
    """
    Compute the Rg from single trajfile.
    """
    
    t = Trajectory.LoadFromLHDF(trajfile)
    Rg = []
    for i in range(len(t['XYZList'])):
        XYZ = t['XYZList'][i,:,:]
        mu = XYZ.mean(0)
        mu = mu.tolist()
        XYZ2 = XYZ-np.tile(mu,(len(XYZ),1))
        Rg.append(((XYZList**2).sum() / n_atoms) ** (0.5))
    return Rg

def checkoutput(Output):
    if Output == None:
	curdir = os.getcwd()
	Output = "%s/Rgs.dat"%curdir
    if os.path.exists(Output) :
	print "%s exists!"%Output
	print "Quit."
	sys.exit()
    return Output

def CalculateRg(Trajfiles,Output='Rgs.dat',returnRgs = False):
    """
    Calculate Radius of gyration for Trajectories.
    
    Trajfiles: can be regular trajectory files ,Gen.lh5 and so on.
    If Trajfiles is None, script will calculate Rgs from trajectories corresponding to ProjectInfo file.
    Output: output file (XXX.dat). 
    The Output default will be set in the scripts and it is './Rgs.dat'.
    """
   
    Output = checkoutput(Output)
  
    Rgs = computeRg(Trajfiles)
	
    print "Save data to %s"%Output
    savetxt(Output,Rgs)
    print "Done."
    
    if returnRgs:
	return Rgs

def CalculateProjectRg(ProjectInfo,Output,returnRgs = False):
    """
    Calculate Radius of gyration for the Project ie. all the Trajectories.
    ProjectInfo: ProjectInfo.h5 file.
    Output: output file (XXX.dat). 
    The Output default will be set in the scripts and it is './Rgs.dat'.
    """
    Output = checkoutput(Output)
   
    if not isinstance(ProjectInfo,str):
	print "Please input the Path to ProjectInfo.h5"
	raise IOError
    print 'Calculating the Rg for each trajectory......' 
    ProjectInfoPath = '/'.join(os.path.realpath(ProjectInfo).split('/')[:-1])
    os.chdir(ProjectInfoPath)
    Trajfiles = []
    ProjectInfo = Serializer.LoadFromHDF(ProjectInfo)
    for i in range(ProjectInfo['NumTrajs']): 
	Trajfiles.append(ProjectInfo['TrajFilePath']+ProjectInfo['TrajFileBaseName']+'%d'%i+ProjectInfo['TrajFileType'])
    Rgs = computeRg(Trajfiles)
    
    print "Save data to %s"%Output
    savetxt(Output,Rgs)
    print "Done."
    if returnRgs:
	return Rgs

def computeRg(trajfiles):
    
    Rgs = []
    i = 1
    try:
	trajfiles = trajfiles.split(' ')
	trajfiles.remove('')
    except (AttributeError,ValueError):
	pass
    for trajfile in trajfiles:
	if os.path.exists(trajfile):
	    print '%d in %d Trajectories'%(i,len(trajfiles)),trajfile
	    Rg = Reference_Rg(trajfile)
	    Rgs.append(Rg)
	    i += 1
	else:
	    print "Can't find trj file:%s",trajfile
	    sys.exit()
    NumTrajs = len(Rgs)
    LongestTraj = max([len(Rgs[i]) for i in range(NumTrajs)])
    Rgs_array = -1*np.ones((NumTrajs,LongestTraj))
    for i in range(NumTrajs):
	Rgs_array[i][:len(Rgs[i])] = Rgs[i][:]
	
    return Rgs_array

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-f','--Trajlikefile',help=''' Path of a Trajectory like file(.lh5).''')
    parser.add_argument('-o','--Output', help='''Output file. default = Rgs.dat''',default='Rgs.dat')
    args = parser.parse_args()
   
    CalculateRg(Trajfiles=args.Trajlikefile,Output=args.Output,returnRgs = False)
    
    