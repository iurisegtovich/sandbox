import os,sys
import string

if len(sys.argv)<2:
    print "Usage: python extractfinalconf.py logfile"
    sys.exit(1)
logfile = sys.argv[1]
fn = open(logfile,'r')
lines = fn.readlines()
fn.close()
goodlines = []

#Old method,from top to bottom
#while 'N-N=' not in lines[0]:
#    lines.pop(0)
#while '\\\\@' not in lines[0]:
#    goodlines.append(lines.pop(0))
#goodlines.append(lines.pop(0))
#print goodlines

#New method, from bottom to top
while '@' not in lines[-1]:
    lines.pop(-1)
while 'N-N=' not in lines[-1]:
    goodlines.append(lines.pop(-1))
goodlines.append(lines.pop(-1))
goodlines.reverse()
#print goodlines

alltext = string.joinfields(goodlines).replace(' ','')
alltext = alltext.replace('\n','')
alltext = alltext.split('\\')
#print alltext

xyz = logfile.replace('log','xyz')
xyzfile = open(xyz,'w')
for line in alltext:
    if line.count(',')==3:
        line = line.replace('0,','')
        line = line.replace(',',' ')
        line = line + '\n'
        xyzfile.writelines(line)

