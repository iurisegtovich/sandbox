import os,sys

for pid in range(6383,6391):
    
    #cmd ="scp tud51931@vav8.chem.temple.edu:/Users/tud51931/projects/FsPeptide_FAH/analysis/first_tica_model/msm/tIC8-lagtime50-states-1200-combined/Data-%d-macro40/tProb.mtx ./Fs-%d-marco40-tProb.mtx"%(pid,pid)
    cmd = "mv Fs-%d-marco40-tProb.mtx Fs-%d-macro40-tProb.mtx"%(pid,pid)
    print cmd
    os.system(cmd)
