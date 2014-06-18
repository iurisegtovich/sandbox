#---------------------------------------
#CHANGE LOG
#---------------------------------------
#Guangfeng Zhou
#Dr.Voelz Lab
#Chemistry Department
#Temple University
#November,1,2012
#---------------------------------------
import os,sys
from scipy import loadtxt

def CalculateRMSError(Populations,Populations_raw):
     
     p = loadtxt(Populations)
     p_raw = loadtxt(Populations_raw)
     if len(p) != len(p_raw):
          print "Populations.dat doesn't match Populations_raw.dat"
          sys.exit()
     else:
          rmserror = (sum([(p[i]-p_raw[i])**2 for i in range(len(p))])/len(p))**0.5
     return rmserror
     
     



           
        
