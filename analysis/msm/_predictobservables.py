#------------------------------------------------
#
#------------------------------------------------
#CHANGE LOG:
#---------------------
#Written by Guangfeng Zhou
#Dr.Voelz Lab
#Chemistry Department
#Temple University
#October,29,2012

#TODO: add some doc strings to help understand the scripts

import os,sys
import numpy as np
from msmbuilder import Serializer,Trajectory
from scipy import savetxt,loadtxt
from scipy.io import mmread
import matplotlib.pyplot as plt
sys.path.append("/Users/tud51931/scripts/gfzhou/")
import HelixCoilTools as hct


__metaclass__ = type

class PredictObservables:
    
    def __init__(self,information,projectfile,populationfile,assignmentfile_fixed,tmatrixfile,rawdatafile):
        
        try:
            self.Info = information
            self.ProjectInfo = Serializer.LoadFromHDF(projectfile)
            self.Population = loadtxt(populationfile)
            self.Assignments = Serializer.LoadFromHDF(assignmentfile_fixed)
            self.Tmatrix = mmread(tmatrixfile)        
            self.StateAssignment = hct.get_StatesAssignments(self.Assignments)
            self.getrawdata(rawdatafile)
        except:
            print "Having trouble with getting required files"
            raise 
        
    def getrawdata(self,rawdatafile):
        
        try:
            d = Serializer.LoadFromHDF(rawdatafile)
            self.RawData = d['Data']
        except:
            try:
                self.RawData = loadtxt(rawdatafile)
            except:
                print "Can not load {}".format(rawdatafile)
                raise 
        
    def predictpopulation(self):
    
        P0 = np.zeros(len(self.Population))
        for data in self.Assignments['Data']:
            P0[data[0]] += 1
        P0 = P0/P0.sum()
        self.PopulationPredicted = []
        for k in range(140):  # tau = 50, so 140*50 = 7000, TODO: change it to a user defined with a default
            self.PopulationPredicted.append(P0)
            P0 *= self.Tmatrix
            
    def PredictData(self):
        
        self.predictpopulation()
        self.calculate_observable_meanvalue_state() # get the mean value of observables for each state    
        self.calculate_observable_meanvalue_step()
        print "Reach here:PredictData"
        self.ObservablePredicted_Step = np.dot(np.array(self.PopulationPredicted), np.array(self.ObservableMeanValue_State).reshape(-1,1))
        self.ObservablePredicted_Step = self.ObservablePredicted_Step.reshape(1,-1)[0]
    
    def calculate_observable_meanvalue_state(self):

        self.Observable_State = {}
        for state in self.StateAssignment.keys():
            for trajid in self.StateAssignment[state].keys():
                for frame in self.StateAssignment[state][trajid]:
                    self.Observable_State.setdefault(state,[]).append(self.RawData[int(trajid)][int(frame)])
                    
        states = [int(i) for i in self.Observable_State.keys()]
        states.sort()
        self.ObservableMeanValue_State = []
        self.ObservableStd_State = []
        for state in states:
            self.ObservableMeanValue_State.append(np.mean(self.Observable_State['%d'%state]))
            self.ObservableStd_State.append(np.std(self.Observable_State['%d'%state]))        
    
    def calculate_observable_meanvalue_step(self):
        
        Counts_ma = np.ma.array(self.RawData,mask=[self.RawData==-1])
        self.ObservableMeanValue_Step = Counts_ma.mean(0)
        self.ObservableStd_State = Counts_ma.std(0)
        
    def PlotPrediction(self,save,output=''):
        
        name = self.Info['ObservableName']
        metric = self.Info['ClusterMetric']
        cutoff = self.Info['Cutoff']
        lagtime = self.Info['Lagtime']
        
        plt.figure()
        plt.plot(np.arange(0,7000,50),self.ObservablePredicted_Step,'ro')# tau = 50, so 140*50 = 7000
        plt.hold(True)    
        plt.plot(range(len(self.ObservableMeanValue_Step)),self.ObservableMeanValue_Step,'b') 
        
        plt.title('%s-steps'%name)
        plt.xlabel('Steps')
        plt.ylabel('%s'%name)
        plt.legend(('%s_Predicted'%name,'%s_RawData'%name),loc = 'upper right')
        if save:
            if output =='':
                figname = '%s_prediction_%sCluster%0.1f_tau%d.png'%(name,metric,cutoff,lagtime)
                output = figname
            plt.savefig(output)
            print "Save to %s"%output        