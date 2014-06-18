import os,sys
from scipy import savetxt,loadtxt

energies = []
states = []
#conffile = '../gaussian_runs/angles_index.dat'
#conf = loadtxt(conffile)
#energy_confs = []
for i in range(100):
    logfile = 'outputs/%d.log'%i
    if os.path.exists(logfile):
        #print conf[i+1]
        fn = open(logfile,'r')
        lines = fn.readlines()
        while 'E(RHF)' not in lines[-1]:
            lines.pop(-1)
        energies.append(float(lines[-1].split()[4]))
    else:
        print "Warning: %s is missing."%logfile
energyfile='./cineromycinB_QM_energies.dat'
#energy_conffile='./Nae-ace-ndm_energyconfs.dat'
print "Wrote: %s"%energyfile
savetxt(energyfile,energies)

