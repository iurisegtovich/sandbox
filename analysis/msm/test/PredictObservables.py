#------------------------------------------------
#
#------------------------------------------------
#CHANGE LOG:
#---------------------
#Written by Guangfeng Zhou
#Dr.Voelz Lab
#Chemistry Department
#Temple University
#October,30,2012

#TODO: add some doc strings to help understand the scripts

import os,sys
import argparse


sys.path.append('Users/tud51931/scripts/gfzhou/test')
from _predictobservables import PredictObservables


parser = argparse.ArgumentParser()
parser.add_argument('-p','--Project',
                    help="""Path to ProjectInfo.h5,default=ProjectInfo.h5""",
                    default="ProjectInfo.h5")
parser.add_argument('-P','--Population',
                    help="""Path to Populations.dat,default=Data/Populations.dat""",
                    default="Data/Populations.dat")
parser.add_argument('-a','--Assignment',
                    help="""Path to Assignments.Fixed.h5,default=Data/Assignments.Fixed.h5""",
                    default="Data/Assignments.Fixed.h5")
parser.add_argument('-t','--Tmatrix',
                    help="""Path to tProb.mtx,default=Data/tProb.mtx""",
                    default="Data/tProb.mtx")
parser.add_argument('-f','--RawData',
                    help="""Path to Raw Data (.dat,.h5)""")
info = {'ObservableName':'Rg','ClusterMetric':'RMSD','Cutoff':4.2,'Lagtime':50}
#parser.add_argument('-i','--Info',
#                    help="""Model Information: Observable name, ClusterMetric, Cutoff, Lagtime. Example:default={'ObservableName':'Nv','ClusterMetric':'RMSD','Cutoff':4.2,'Lagtime':50}""",
#                    default={'ObservableName':'Rg','ClusterMetric':'RMSD','Cutoff':4.2,'Lagtime':50})
parser.add_argument('-n','--name',help="Observable Name")
parser.add_argument('-m','--metric',help="ClusterMetric,default = RMSD",default='RMSD')
parser.add_argument('-c','--cutoff',help="Cutoff,default=4.2",default=4.2)
parser.add_argument('-l','--lagtime',help="Lagtime,default=50",default=50)
args = parser.parse_args()
info['ObservableName'] = args.name
info['ClusterMetric'] = args.metric
info['Cutoff'] = float(args.cutoff)
info['Lagtime'] = int(args.lagtime)

print info,args.Project,args.Population,args.Assignment,args.Tmatrix,args.RawData

Prediction=PredictObservables(information=info,
                              projectfile=args.Project,
                              populationfile=args.Population,
                              assignmentfile_fixed=args.Assignment,
                              tmatrixfile=args.Tmatrix,
                              rawdatafile=args.RawData)
Prediction.PredictData()
Prediction.PlotPrediction(save=True)
