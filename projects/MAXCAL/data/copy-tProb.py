import os,sys
for i in range(6383,6391):
    #cmd = "scp tud51931@vav8.chem.temple.edu:/Users/tud51931/projects/FsPeptide_FAH/analysis/first_tica_model/msm/tIC8-lagtime50-states-1200-combined/Data-%i-macro40/tProb.npy ./tProb-%d.npy"%(i,i)
    #cmd = "scp tud51931@vav8.chem.temple.edu:/Users/tud51931/projects/FsPeptide_FAH/analysis/first_tica_model/msm/tIC8-lagtime50-states-1200-combined/Data-%i-macro40/tCounts.npy ./tCounts-%d.npy"%(i,i)
    cmd = "scp tud51931@vav8.chem.temple.edu:/Users/tud51931/projects/FsPeptide_FAH/analysis/first_tica_model/msm/tIC8-lagtime50-states-1200-combined/Data-%i-macro40/Nhs_state.dat ./Nhs_state-%d.dat"%(i,i)

    os.system(cmd)
