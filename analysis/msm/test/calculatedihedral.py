from msmbuilder import Trajectory,Serializer
from msmbuilder.geometry import dihedral
from scipy import savetxt
import numpy as np
import matplotlib.pyplot as plt


def is_helical(phiangle, psiangle):
    """Returns True if the phi/psi angles(degrees) respresent
    a helical conformation."""

    phirange = [-60.0-30.0, -60.0+30.0]
    psirange = [-47.0-30.0, -47.0+30.0]

    if phiangle > phirange[0] and phiangle < phirange[1]:
        if psiangle > psirange[0] and psiangle < psirange[1]:
            return True
    return False

def get_StatesAssignments(AssignmentFiles):
    """
    StatesAssignments {'state':{'Trajectory':[Frames]}}
    """
    A = AssignmentFiles
    S = {}
    try:
        for trajid,data in zip(A['TrajID'],A['Data']):
            for i in len(data):
                if data[i] == -1:
                    continue
                else:
                    S.setdefault('%d'%data[i],{}).setdefault('%d'%trajid,[]).append(i)
        return S
    except KeyError:
        # Assignments file which doesn't has key 'TrajID' should be regular Assignments file instead of new assignments files created for bootstrap. Then modify the assignments file i.e, create key 'TrajID'.
        A['TrajID'] = list(range(len(A['Data'])))
        get_StatesAssignments(A)
        
def ComputeDihedralsFromAssignmentFile(AssignmentFiles):
    """
    Compute dihedrals from Assignment files(Assignment_fixed) 
    This is used for compute dihedrals for mapped states.
    
    SA is StatesAssignments
    SA_dihedral
    """
    SA = get_StatesAssignments(AssignmentsFiles)
    states = [SA.keys()]
    trajectory = {}
    for state in states:
        trajids = [SA['%d'%state].keys()]
        for trajid in trajids:
            trajectory['XYZList'] = get_Trajectory_frame(trajid,SA['%d'%state]['%d'%trajid])            
            dihedrals = ComputeDihedralsFromTrajectory(trajectory)
            
            
def ComputeDihedralsFromTrajectory(Trajectory):
    """
    Compute dihedrals directly from trajectory files.
    """
    
    traj = {}
    traj['XYZList'] = Trajectory['XYZList']/1.0 #convert XYZList to float type
    indices = dihedral.get_indices(Trajectory,'phi/psi')
    dihedrals = dihedral.compute_dihedrals(traj,indices,'phi/psi')
    
    return dihedrals

def get_Trajectory_frame(trajid,frames):
    """
    Get trajectory frames.
    From Trajectory file(traj_.lh5) get the frame.
    """
    Path = "/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/sourcedata/Trajectories"
    traj = Serializer.LoadFromHDF('%s/trj%d.lh5'%(Path,trajid))
    
    return [traj['XYZList'][i] for i in frames]
    
def CreateTrajFileWithDihedralList(Trajectory):
    """
    Create Trajectory file with Dihedral list
    """
    T = Trajectory
    dihedrals = ComputeDihedralsFromTrajectory(T)
    T['Dihedrals'] = dihedrals
    names = T['SerializerFilename'].split('.')
    names[0] += '_Dih'
    T['SerializerFilename'] = names[0]+'.'+'h5'
    print "Save to %s"%T['SerializerFilename']
    T.SaveToHDF(T['SerializerFilename'])

def CreateTrajFileWithHCstrings(Trajectory):
    """
    Create Trajectory file with HC strings
    """    
    T = Trajectory
    HCs = []
    dihedrals = ComputeDihedralsFromTrajectory(T)
    for i in range(len(dihedrals)):
        hc = ''
        for j in range(len(dihedrals[i])/2):
            if is_helical(dihedrals[i][j],dihedrals[i][j+len(dihedrals[i])/2]):
                hc +='h'
            else:
                hc +='c'
        HCs.append(hc)
    T['HCs'] = HCs
    names = T['SerializerFilename'].split('.')
    names[0] += '_hc'
    T['SerializerFilename'] = names[0]+'.'+'h5'
    print "Save to %s"%T['SerializerFilename']
    T.SaveToHDF(T['SerializerFilename'])    

def count_Helix(Trajectory):
    """
    Read in 'HCs' from modified trajectory files(trj_hc.h5) 
    """
    T = Trajectory
    try:
        Hcount = [frame.count('h') for frame in T['HCs']]
        return Hcount
    except KeyError:
        print "Can't find key 'HCs', please use trj_hc.h5"
        raise KeyError
            

ProjectInfo = Serializer.LoadFromHDF('/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/sourcedata/ProjectInfo.h5')
Counts = -1*np.ones((ProjectInfo['NumTrajs'],max(ProjectInfo['TrajLengths'])))
print Counts.shape

for i in range(0,93):
    T = Trajectory.LoadFromHDF('/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/sourcedata/Trajectories/trj%d_hc.h5'%i)
    Hcount = count_Helix(T)
    Counts[i,:len(Hcount)] = Hcount[:]

Counts_ma = np.ma.array(Counts,mask=[Counts==-1])
H_mean = Counts_ma.mean(0)
H_std = Counts_ma.std(0)
print H_mean

plt.figure()
plt.plot(range(len(H_mean)),H_mean)
plt.savefig('numhelix-rawdata.png')
    

