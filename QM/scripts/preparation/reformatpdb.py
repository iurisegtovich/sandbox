import os,sys

for i in range(100):
    pdb = open('Gen%d_minimized.pdb'%i,'r')
    lines = pdb.readlines()
    pdb.close()
    newpdb = open('%d.pdb'%i,'w')
    newlines = lines[2:len(lines)-2]
    newpdb.writelines(newlines)
    newpdb.writelines('END')
    newpdb.close()
