#A short demonstration to show how to use buildmsm.py
#Guangfeng Zhou
#Voelz Lab
#Chemistry Department
#Temple University

import os,sys

#Use -t, you can specify one specific directory contains all the data you need to build MSM.
os.system("python buildmsm.py -t /Users/tud51931/projects/MSM/msm/ff03ERR-hybridkcenter/RMSDCluster3.0 -c 3 -l 50 100 150 -m 'rmsd hybrid'")
#Use -i, you can specify a parent directory contains subdirectories then for each subdirectories, you can build them under different cutoffs -- use -c cutoffs. But this requires certain folder structure.
os.system("python buildmsm.py -i /Users/tud51931/projects/MSM/msm/ff03ERR-hybridkcenter/ -c 3 3.5 4 4.5 5 -l 50 100 150 -m 'rmsd hybrid'")