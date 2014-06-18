from math import ceil
from scipy.io import mmread


finalmatrix = '/Users/tud51931/projects/MSM/msm/ff03ERR-hybridkcenter/RMSDCluster4.0/Dataforsurprisal/tCounts.mtx'
fix = mmread(finalmatrix)
fix = ceil(fix*100.0/fix.sum())
print fix.todense()
