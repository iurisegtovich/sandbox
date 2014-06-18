import os,sys

for i in range(100):
    pdbfile = '/Users/tud51931/projects/QM/albocycline/albo-cineromycinB/minimized_pdbs_for_QM/output/%d.pdb'%i
    fn = open(pdbfile,'r')
    lines = fn.readlines()
    fn.close()
    for j in range(len(lines)):
        lines[j] = lines[j].replace('HETATM','ATOM  ')
    newpdbfile = pdbfile.replace('.pdb','.fixed.pdb')
    newfn = open(newpdbfile,'w')
    newfn.writelines(lines)
    print "Wrote: %s"%newpdbfile
    newfn.close()

