import os,sys
from scipy import savetxt,loadtxt

energies = []
states = []
#conffile = '../gaussian_runs/angles_index.dat'
#conf = loadtxt(conffile)
#energy_confs = []
for i in range(100):
    logfile = 'output/%d.log'%i
    if os.path.exists(logfile):
        #print conf[i+1]
        fn = open(logfile,'r')
        lines = fn.readlines()
        for line in lines:
            if 'E(RB3LYP)' in line:
                #energy_confs.append(conf[i])
                states.append(i)
                energies.append(float(line.split()[4]))
energyfile='./cineromycinB_QMenergies.dat'
#energy_conffile='./Nae-ace-ndm_energyconfs.dat'
print "Wrote: %s"%energyfile
savetxt(energyfile,energies)

