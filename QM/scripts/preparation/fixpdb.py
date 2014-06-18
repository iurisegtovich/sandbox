import os,sys

for i in range(100):
    pdbfile = '%d.pdb'%i
    fn = open(pdbfile,'r')
    lines = fn.readlines()
    fn.close()
    for i in range(len(lines)):
        lines[i] = lines[i].replace('HETATM','ATOM  ')
    newpdbfile = pdbfile.replace('.pdb','.fixed.pdb')
    newfn = open(newpdbfile,'w')
    newfn.writelines(lines)
    print "Wrote: %s"%newpdbfile
    newfn.close()

