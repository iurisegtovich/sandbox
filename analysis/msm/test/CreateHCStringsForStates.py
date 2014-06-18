import sys,os
sys.path.append('/Users/tud51931/scripts/gfzhou')
from msmhcanalysis import create_hcstrings_states

cutoff = 5.2
metrics = 'dihedral'

if metrics.lower() == 'dihedral':
    Path = "/Users/tud51931/projects/MSM/msm/ff03-dihedralhybrid/"
    metrics = 'Dihedral'
elif metrics.lower() == 'rmsd':
    Path = "/Users/tud51931/projects/MSM/msm/ff03-hybridkcenter/"
    metrics = 'RMSD'
Path = os.path.join(Path,'%sCluster%0.1f'%(metrics,cutoff))
AssignmentFile = os.path.join(Path,"Data","Assignments.h5")


create_hcstrings_states(Assignments=AssignmentFile,outfile = 'HCstrings_states_%s%0.1f.txt'%(metrics,cutoff))