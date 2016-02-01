#SurprisalBasedSampling.py

import os,sys
sys.path.append('/Users/tud51931/voelzlab/surprisal/')
import scipy
import random
import commands
import argparse
import numpy as np
from scipy import savetxt,loadtxt
from scipy.io import mmread,mmwrite
from Surprisal import SurprisalCalculator,get_populations
from msmbuilder import MSMLib
from msmbuilder.io import saveh,loadh
import time





def _get_most_jsd_state(allstates_jsd):
    
    StateID = np.where(allstates_jsd==max(allstates_jsd))[0][0]
    return StateID

def _get_most_jsd_var_state(allstates_jsd_var):

    StateID = np.where(allstates_jsd_var==max(allstates_jsd_var))[0][0]
    return StateID    

def _compute_simulation_time_from_sparse_matrix(sparse,lagtime,n_trajs,frame_interval=100.0):
    
    time = round((sparse.sum()+lagtime*n_trajs)*frame_interval/10**6,6) #micro seconds
    print "simulation time",time
    return time

class SurprisalBasedAdaptiveSampling:
    
    def __init__(self,sparse1,sparse2,criterion,lagtime,n_trajs,n_frames):
        
        self.sparse1 = sparse1
        self.sparse2 = sparse2
        self.criterion = criterion
        self.lagtime = lagtime
        self.n_trajs = n_trajs
        self.n_frames = n_frames
        self.criterion_counter = 0  
        self.populations1 = get_populations(self.sparse1, 'mle')[0] 
        print self.populations1.shape
        self.populations2 = []        
        
    def check_criterion(self):
        if isinstance(self.criterion,int):
            if self.criterion >= self.criterion_counter:
                self.criterion_counter += 1
                return True
            else:
                return False
        else:
            print "Not supported criterion type."
            raise TypeError
        
    def _fix_sparse(self,totalmatrix=None):
        """
        Add 1 count for each transition
        
        size = range(self.sparse2.shape[0])
        x,y = np.meshgrid(size,size)
        data = np.ones(len(x.flatten()),dtype=int)
        fix = scipy.sparse.coo_matrix((data,(x.flatten(),y.flatten())),shape=self.sparse2.shape)
        self.sparse2 = self.sparse2 +fix
        """
        """        
        #Add 1 count to self transition
        size = range(self.sparse2.shape[0])
        data = np.ones(len(size),dtype=int)
        fix = scipy.sparse.coo_matrix((data,(size,size)),shape=self.sparse2.shape)
        print fix.todense()
        self.sparse2 = self.sparse2 +fix
        """
        """        
        # Use fractional total count matrix
        if totalmatrix != None:
            fix = mmread(totalalmatrix)
            self.sparse2 = self.sparse2 + fix*100.0/fix.sum()
        """
        # ceil franctional counts to 1
        if totalmatrix != None:
            fix = mmread(totalmatrix)
            fix = fix/fix.sum()
            fix = fix.asformat('csr')
            for i in range(fix.shape[0]):
                for j in range(fix.shape[0]):
                    if fix[i,j] != 0:
                        fix[i,j] = np.ceil(fix[i,j])
            fix = fix.asformat('coo')
            self.sparse2 = self.sparse2 + fix        
        
        
    def pseudosampling(self, countslib ,algorithm):
        
        simulation_time=[]
        jsds = []
        jsd_vars = []
        all_sparse=[]
        stateid = -1
        if os.path.exists(countslib):
            mtxfiles = os.path.join(countslib,'*.mtx')
            all_counts_files = commands.getoutput('ls %s'%mtxfiles).split('\n')
            for i,fn in enumerate(all_counts_files):
                sparse = mmread(fn)
                all_sparse.append(sparse)
                print "Loading in %d/%d matrices."%(i+1,len(all_counts_files))        
        while True:
            
            if (self.criterion_counter % 10 == 0):
                print "Using MLE"
                if self.criterion_counter == 0:
                    print "Fixing"
                    self._fix_sparse(totalmatrix='/Users/tud51931/projects/MSM/msm/ff03ERR-hybridkcenter/RMSDCluster4.0/Dataforsurprisal/tCounts.mtx')
                self._compute_allstates_jsd(popmethod='retain_first_mle')
                
            else:
                print "Using retain"
                self._compute_allstates_jsd(popmethod='retain')                
            
            print "start jsd_var"
            self._compute_allstates_jsd_var(stateid,n_bootstraps=50)       
            simulation_time.append(_compute_simulation_time_from_sparse_matrix(self.sparse2,self.lagtime,n_trajs=self.n_trajs,frame_interval=100.0))
            print "finish jsd_var"
            jsds.append(np.sum(self.allstates_jsd))
            jsd_vars.append(np.sum(self.allstates_jsd_var))
            print "step %d of %d"%(self.criterion_counter,self.criterion)
            
            if not self.check_criterion(): # check whether continue or not.
                break
            if algorithm.lower() == 'alg1':
                stateid = _get_most_jsd_state(self.allstates_jsd)
                print stateid
            elif algorithm.lower() == 'alg2':
                stateid = _get_most_jsd_var_state(self.allstates_jsd_var)
                print stateid
            
            print "JSD:",np.sum(self.allstates_jsd)
            print "JSD_var:",np.sum(self.allstates_jsd_var)
            
            pseudosparse = self.generate_pseudo_sparse(self.n_trajs,stateid,all_sparse,self.n_frames)
            self.sparse2 = self.sparse2 + pseudosparse
        
        return simulation_time,jsds,jsd_vars
    
    
    def generate_pseudo_sparse(self,numberoftrajs,stateid,all_sparse,frames):
        
        initialize=False
        for traj in range(numberoftrajs):
            sparse= random.choice(all_sparse)
            counts = sparse.toarray()[stateid]
            if counts.sum() == 0:
                prob = np.zeros(len(counts))
            else:
                prob = counts/counts.sum()
            if not initialize:
                pseudo_counts = np.zeros(len(counts))
                initialize = True
            pseudo_counts += np.random.multinomial(frames/self.lagtime,prob)
        row = [stateid]*len(counts)
        col = range(len(counts))
        pseudo_sparse = scipy.sparse.coo_matrix((pseudo_counts,(row,col)),shape=sparse.shape)
        return pseudo_sparse
    
    def _compute_allstates_jsd(self,popmethod = 'counts'):
       
        if popmethod == 'retain_first_mle': 
            self.obj = SurprisalCalculator(self.sparse1, self.sparse2, popmethod, self.populations1)      
            self.populations2 = self.obj.populations2
            self.allstates_jsd = self.obj.calculate_all_jsd()

        if popmethod == 'retain':
            self.obj = SurprisalCalculator(self.sparse1, self.sparse2, popmethod, self.populations1, self.populations2)
            self.allstates_jsd = self.obj.calculate_all_jsd()
        
    
    def _compute_allstates_jsd_var(self,stateid,n_bootstraps=50):
        
        if (stateid ==-1) :
            print "all_jsd_var"
            self.allstates_jsd_var = self.obj.calculate_all_jsd_var()
            print "Finish all_jsd_var"
        else:
            self.allstates_jsd_var[stateid] = self.obj.calculate_jsd_var(stateid, n_bootstraps=50)


if __name__ == '__main__':
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-oc','--originalcounts',help="Path to Original tCounts.mtx file",required=True)
    parser.add_argument('-a','--assignments',help="Path to Assignments.h5 file",required=True)
    parser.add_argument('-f','--frames',help="The number of frames in each of the new trajecotries.Default = 50 (5ns)",default=50,type=int)
    parser.add_argument('-o','--output',help="Pseudo Trajectory Assignment file.",required=True,type=str)
    parser.add_argument('-v','--verbose',help="increase output verbosity",action="store_true")
    args = parser.parse_args()
    sparse1 = mmread(args.originalcounts)
    """

    sparse1 = mmread('/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster4.0/Dataforsurprisal/tCounts.mtx')
    sparse2 = mmread('/Users/tud51931/projects/MSM/msm/ff03ERR-hybridkcenter/RMSDCluster4.0/Dataforsurprisal/SlicedAssignments/UnMappedCounts/tCounts-100.mtx')
    countslib = '/Users/tud51931/projects/MSM/msm/ff03ERR-hybridkcenter/RMSDCluster4.0/Dataforsurprisal/100trajs/bs0/UnMappedCounts'
    project = SurprisalBasedAdaptiveSampling(sparse1,sparse2,criterion=2000,lagtime=50,n_trajs=93,n_frames=50)
   
    simutime,jsds,jsd_vars = project.pseudosampling(countslib,algorithm='alg2')
    data = np.zeros((3,len(simutime)))
    data[0],data[1],data[2] = simutime,jsds,jsd_vars
    print simutime,jsds,jsd_vars
    savetxt('ERRdata_JSD_PseudoSampling_ALG2_MLE.dat',data)

