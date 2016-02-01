import os,sys
for i in range(100,800,100):
    cmd = 'python CalculateJSD_BootStrapped.py -oc /Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/RMSDCluster4.0/Dataforsurprisal/tCounts.mtx -bc /Users/tud51931/projects/MSM/msm/ff03ERR-hybridkcenter/RMSDCluster4.0/Dataforsurprisal/%dtrajs/bs0/UnMappedCounts -m retain_first_mle -f -o JSDs_ERR_RMSD4.0_MLE_%d.trajs.dat '%(i,i)
    print cmd
    os.system(cmd)
