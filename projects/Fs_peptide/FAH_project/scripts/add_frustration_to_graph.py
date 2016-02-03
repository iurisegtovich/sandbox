import sys, os, glob
import numpy as np
from scipy import loadtxt

# read in the DOT graph
fin = open('Graph.dot','r')
lines = fin.readlines()
fin.close()

# read in the frustrations
#gin = open('../../../gb1-fah/figXX-MFPTs/frustration_wt.dat','r')
gin = open('../../../gb1-fah/figXX-MFPTs/frustration_wt-4x.dat','r')
glines = gin.readlines()[1:]  # pop the header
gin.close()

# read in the populations
pops = loadtxt('Populations.dat')

# read in the JSD values
#jsds = np.zeros(150)
#jsd_data = loadtxt('sorted_jsd_wt_tz4.dat')
#for i in range(jsd_data.shape[0]) :
#    state = jsd_data[i,0]
#    jsds[state] = jsd_data[i,1]


frustrations = np.zeros(150)

for line in glines:
    state = int(line.split()[0])
    f = float(line.split()[1])
    frustrations[state] = f 

newlines = []
for line in lines:

    s = line.strip()
    if (s[-2:] == '];') and (s.count('->') == 0):
       state = int(s.split(' ')[0])
       #snew = s.replace('];', ', frust="%f", population="%f", jsd="%f", logjsd="%f"];'%(frustrations[state], pops[state], jsds[state], np.log(jsds[state])) )
       snew = s.replace('];', ', frust="%f", population="%f"];'%(frustrations[state], pops[state]) )
       newlines.append( snew+'\n' )

    else:
       newlines.append( line )

# Write output file
outfile = 'Graph_frust_wt-4x.dot'
fin = open(outfile,'w')
fin.writelines(newlines)
fin.close()

for line in newlines:
    print line.strip()

print
print 'Wrote', outfile







