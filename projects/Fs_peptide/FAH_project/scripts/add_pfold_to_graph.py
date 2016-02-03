import sys, os, glob
import numpy as np
from scipy import loadtxt

# read in the DOT graph

for p_id in range(6383,6391):
    fin = open('Graph-{pid}-macro40-netflux.dot'.format(pid=p_id),'r')
    lines = fin.readlines()
    fin.close()

    # read in the frustrations
    #gin = open('../../../gb1-fah/figXX-MFPTs/frustration_wt.dat','r')
    committors_fn = 'Data-{pid}-macro40/committors.npy'.format(pid=p_id)
    committors = np.load(committors_fn)

    # read in the populations
    pops = loadtxt('Data-{pid}-macro40/Populations.dat'.format(pid=p_id))
 
    freeE = loadtxt('Data-{pid}-macro40/Free_energies.dat'.format(pid=p_id))

    # read in the JSD values
    #jsds = np.zeros(150)
    #jsd_data = loadtxt('sorted_jsd_wt_tz4.dat')
    #for i in range(jsd_data.shape[0]) :
    #    state = jsd_data[i,0]
    #    jsds[state] = jsd_data[i,1]

    newlines = []
    for line in lines:

        s = line.strip()
        if (s[-2:] == '];') and (s.count('->') == 0):
           state = int(s.split(' ')[0])
           #snew = s.replace('];', ', frust="%f", population="%f", jsd="%f", logjsd="%f"];'%(frustrations[state], pops[state], jsds[state], np.log(jsds[state])) )
           snew = s.replace('];', ', pfold="%f", population="%f", free_energies="%f"];'%(committors[state], pops[state], freeE[state]) )
           newlines.append( snew+'\n' )

        else:
           newlines.append( line )

    # Write output file
    outfile = 'Graph-{pid}-macro40-netflux-pfold.dot'.format(pid=p_id)
    fin = open(outfile,'w')
    fin.writelines(newlines)
    fin.close()

    for line in newlines:
        print line.strip()

    print
    print 'Wrote', outfile


