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
    if isinstance(A,str):
        A = Serializer.LoadFromHDF(A)
    S = {}
    try:
        for trajid,data in zip(A['TrajID'],A['Data']):
            for i in range(len(data)):
                if data[i] == -1:
                    continue
                else:
                    S.setdefault('%d'%data[i],{}).setdefault('%d'%trajid,[]).append(i)
        return S
    except KeyError:
        # Assignments file which doesn't has key 'TrajID' should be regular Assignments file instead of new assignments files created for bootstrap. Then modify the assignments file i.e, create key 'TrajID'.
        A['TrajID'] = list(range(len(A['Data'])))
        return get_StatesAssignments(A)

        
#def ComputeDihedralsFromAssignmentFile(AssignmentFiles):
#    """
#    Compute dihedrals from Assignment files(Assignment_fixed) 
#    This is used for compute dihedrals for mapped states.
#    
#    SA is StatesAssignments
#    SA_dihedral
#    """
#    SA = get_StatesAssignments(AssignmentsFiles)
#    states = SA.keys()
#    trajectory = {}
#    for state in states:
#        trajids = SA['%d'%state].keys()
#        for trajid in trajids:
#            trajectory['XYZList'] = get_Trajectory_frame(trajid,SA['%d'%state]['%d'%trajid])            
#            dihedrals = ComputeDihedralsFromTrajectory(trajectory)
            
            
def ComputeDihedralsFromTrajectory(Trajectory):
    """
    Compute dihedrals directly from trajectory files.
    """
    
    traj = {}
    traj['XYZList'] = Trajectory['XYZList'] #convert XYZList to float type
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
    
def ConvertDihedralsToHCStrings(dihedrals):
    HCs = []
    for i in range(len(dihedrals)):
        hc = ''
        for j in range(len(dihedrals[i])/2):
            if is_helical(dihedrals[i][j],dihedrals[i][j+len(dihedrals[i])/2]):
                hc +='h'
            else:
                hc +='c'
        HCs.append(hc)
    return HCs
    
def CreateTrajFileWithHCstrings(Trajectory):
    """
    Create Trajectory file with HC strings
    """    
    T = Trajectory
    dihedrals = ComputeDihedralsFromTrajectory(T)
    HCs = ConvertDihedralsToHCStrings(dihedrals)
    T['HCs'] = HCs
    names = T['SerializerFilename'].split('.')
    names[0] += '_hc'
    T['SerializerFilename'] = names[0]+'.'+'lh5'
    print "Save to %s"%T['SerializerFilename']
    T.SaveToLHDF(T['SerializerFilename'])    

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

def compute_numhelix_trajectory(Trajectory):
    """
    compute number of helix from original traj file i.e. trj%d.lh5    
    """
    dihedrals = ComputeDihedralsFromTrajectory(Trajectory)
    numhelix = []
    for i in range(len(dihedrals)):
        helixcount = 0
        for j in range(len(dihedrals[i])/2):
            if is_helical(dihedrals[i][j],dihedrals[i][j+len(dihedrals[i])/2]):
                helixcount += 1
        numhelix.append(helixcount)
    return numhelix

def compute_numhelix_states(StatesAssignments):
    """
    Compute the average number of helix for all states.
    Need Path to trj_hc.h5
    """
    
    SA = StatesAssignments
    states = SA.keys()
    numhelix_states = {}
    n = 0
    for state in states:
        n +=1
        print "Compute number of helix for state %d/%d"%(n,len(states))
        TrajID = SA[state].keys()
        numhelix_state = []
        for trajid in TrajID:
            T = {}            
            TrajFile = '/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/sourcedata/Trajectories/trj%s_hc.h5'%trajid
            Traj = Serializer.LoadFromHDF(TrajFile)
            T['HCs'] = [Traj['HCs'][i] for i in SA[state][trajid]]
            numhelix_state += count_Helix(T)
        numhelix_states[state] = numhelix_state
    
    return numhelix_states
        
def test():
    assignments = Serializer.LoadFromHDF("/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster4.2/Data/Assignments.h5")
    StatesAsi = get_StatesAssignments(assignments)
    NumHelix_states = compute_numhelix_states(StatesAsi)
    #print "NumHelix_states",NumHelix_states
    #savetxt('NumHelix_states',NumHelix_states)
    states = [int(i) for i in NumHelix_states.keys()]
    states.sort()
    mean_numhelix_states = []
    std_numhelix_states = []
    for state in states:
        mean_numhelix_states.append(np.mean(NumHelix_states['%d'%state]))
        std_numhelix_states.append(np.std(NumHelix_states['%d'%state]))
    
    plt.figure()
    plt.errorbar(states,mean_numhelix_states,std_numhelix_states)
    plt.xlabel("State ID")
    plt.ylabel("Number of Helix")
    plt.savefig("Numhelix_states")
    plt.show()    
    
def test1():
    """
    This test shows how to get the number of helix from a trajectory.
    """
    traj = Trajectory.LoadFromLHDF('/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/sourcedata/Trajectories/trj34.lh5')
    numhelix = compute_numhelix_trajectory(traj)
    print len(numhelix)
    
def test2():
    """
    This test shows how to create new trj files with hc strings.
    """
    path ="/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/sourcedata/Trajectories"
    for i in range(0,100):
        Trajfile = "%s/trj%d.lh5"%(path,i)
        if os.path.exists(Trajfile):
            T = Trajectory.LoadFromLHDF(Trajfile)
            CreateTrajFileWithHCstrings(T)
    print "Done."
            
if __name__ == '__main__' :
    test2()