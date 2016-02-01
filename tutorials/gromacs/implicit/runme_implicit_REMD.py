#!/usr/bin/env python

import sys, os
import random, math, string, tempfile
import numpy as np

usage = """runme_explicit_REMD.py GROFILE TOPFILE NPROC NREPLICAS 

   Prepares a gromacs simulation from a PDBFILE,
   and runs on NPROC number of processors.

   The number of temperatures will be set to the
   number of processors, exponentially-spaced in
   the range 300-600K"""

if len(sys.argv) < 3:
    print usage
    sys.exit(1)

GroFile = "conf.gro"
TopFile = "equilibrated.top"
nproc = 12 # The number of processors - must be a multiple of nreplicas
nreplicas = 12 # The number of replicas
ntemps = nreplicas

UseOwlsnest = True #False
SkipEquil = False  #True

# DEBUG flag:
DEBUG = False #True
TESTING = False  # if True, print cmds but dpn't run them

def cleanup():
    """Delete files not needed."""

    os.system('rm *.gro; rm *.ndx; rm *.log')
    os.system('rm *.top; rm *.tpr')
    os.system('rm step*')
    os.system('rm ./#*')

def run_cmd(cmd, testing=False):
    """Execute command-line command."""

    print '>>', cmd
    if not testing:
        os.system(cmd)

def write_production_mdp(prefix, rep, temp):
    """Write a production mdp file "prefix"_<rep>.mdp with temperature temp."""

    fn = prefix + '_%d.mdp'%rep
    fout = open(fn, 'w')
    content = """title                    = prod_GBSA
cpp                      = /lib/cpp
include                  = -I../top
define                   =
integrator               = sd
dt                       = 0.002
nsteps                   = 50000000   ; 100 ns
nstxout                  = 50000
nstvout                  = 50000
nstlog                   = 50000
nstxtcout                = 50000

implicit_solvent = GBSA
gb_algorithm     = OBC
nstgbradii       = 1
gb_epsilon_solvent = 80.0
sa_algorithm    = Ace-approximation

comm_mode = ANGULAR

rgbradii         = 0.9
coulombtype     = Cut-off
rvdw            = 0.9
rlist           = 0.9
rcoulomb        = 0.9

; CONSTRAINTS
constraints     = hbonds
constraint_algorithm = LINCS

; other stuff
bd_fric        = 1.0
nstcomm                  = 10
comm_grps                = System
xtc_grps                 = Protein
energygrps               = System
nstlist                  = 10
ns_type                  = grid
tc-grps                  = System
tau_t                    = 0.0109

ref_t                    = {temp}
compressibility          = 4.5e-5
ref_p                    = 1.0
gen_vel                  = yes
gen_temp                 = {temp}
gen_seed                 = 1729
pbc                      = no
    """.format(temp=temp)
    fout.write(content)
    fout.close()
    return fn


#-----------------------------------
# MAIN
#-----------------------------------

if UseOwlsnest:

    #nthreads = 8     # devel queue
    #queue = 'devel'
    #walltime = '0:30:00'

    nthreads = 12    # normal queue
    queue = 'normal'
    walltime = '0:30:00'

    nnodes = nproc


else:
    #nthreads = 8     # for vav8 8-core MacPro
    nthreads = 12    # for vav2 12-core MacPro
    #nthreads = 4      # for vav1, vav5,6,7 4-core Mac Mini Server
    #nthreads = 1      # for small systems

# cleanup()  # cleanup from last time


# make exponentially-spaced temps
mintemp = 300.
maxtemp = 600.
dlogtemp = (np.log(maxtemp) - np.log(mintemp))/(ntemps-1)
temps = np.exp(np.arange(np.log(mintemp), np.log(maxtemp)+dlogtemp, dlogtemp ))
temps = np.round(temps)[0:ntemps]
print 'temps', temps

"""
# editconf -- generate and/or redefine the box coordinates 
run_cmd('editconf -f %s -o Conf.gro -bt octahedron -d 1.0'%GroFile, testing=TESTING)
 
#fill the box with water
run_cmd('genbox -cp Conf.gro -cs spc216.gro -o conf_solvated.gro', testing=TESTING) 

if (0):
    #make a new topology
    run_cmd('pdb2gmx -f conf_solvated.gro -o NEW_conf_solvated.gro -p tNEW_conf_solvated.top -ignh', testing=TESTING)
else:
    pass  # we're using pep5_solvated.top  which we edited by hand 

# Make an index file:
run_cmd( "echo 'name 1 Peptoid\nq\n' | make_ndx -f conf_solvated.gro -o index.ndx", testing=TESTING)

GroFile = 'conf_solvated.gro'
TopFile = 'pep5_solvated.top'

# Add counterions at ~100 mM to neutralize the system
run_cmd( 'grompp -f mdp/minimize.mdp -c %s -p %s -n index.ndx -o minimize.tpr'%(GroFile, TopFile), testing=TESTING )
run_cmd( "echo '7\n' | genion -s minimize.tpr -n index.ndx -p pep5_solvated.top -o solvent_ions.gro -pname NA -nname CL  -neutral -conc 0.1", testing=TESTING )
run_cmd( "echo 'name 1 Peptoid\nq\n' | make_ndx -f solvent_ions.gro -o index.ndx", testing=TESTING )

# Minimize
run_cmd( 'grompp -f mdp/minimize.mdp -c solvent_ions.gro -p pep5_solvated.top -n index.ndx -o minimize2.tpr', testing=TESTING )
run_cmd( 'mdrun -v -nt 8 -s minimize2.tpr -c solvent_ions_minimized.gro', testing=TESTING )

# Equilibrate with pressure coupling on (Berendsen)
run_cmd( 'grompp -f mdp/equil.mdp -c solvent_ions_minimized.gro -p pep5_solvated.top -n index.ndx -o equil.tpr', testing=TESTING )
run_cmd( 'mdrun -v -nt 8 -s equil.tpr -c solvent_ions_equilibrated.gro', testing=TESTING )

"""

GroFile = 'solvent_ions_equilibrated.gro'
TopFile = 'pep5_solvated.top'

# Set up production runs  for all temps
for rep in range(len(temps)):
    temp = temps[rep]
    mdpfilename = write_production_mdp("prod", rep, temp)
    run_cmd( 'grompp -f %s -c %s -p %s -n index.ndx -o prod_%d.tpr'%(mdpfilename, GroFile, TopFile, rep), testing=TESTING )
    
if UseOwlsnest:
    # write a qsub script to the rundir
    qsubfile = 'qsub.sh'
    ppn = nthreads
    fout = open(qsubfile,'w')
    fout.write("""#!/bin/sh
#PBS -l walltime=%s
#PBS -N peptoid_ex_REMD  
#PBS -q %s 
#PBS -l nodes=%d:ppn=%d
#PBS -o peptoid_ex_REMD 
#PBS 

cd $PBS_O_WORKDIR

module load gromacs
module load openmpi
#time mdrun_mpi -v -s prod_.tpr -multi 12 -replex 5000 -c conf_prod.gro
mpirun mdrun_mpi -v -nt %d -s prod_.tpr -multi %d -replex 5000 -c conf_prod.gro
"""%(walltime, queue, nnodes, ppn, ppn, nnodes)   )
    fout.close()     # note that we don't need "mpirun -np 12" -- the cluster automatically figures this out
    
    #if not TESTING:
    #    os.system('qsub qsub.sh')
    
else:
    # print the mdrun_mpi command we need to run
    print 'Run this:'
    run_cmd( 'mpirun -np %d mdrun_mpi -v -s prod_.tpr -multi %d -replex 500 -c conf_prod.gro'%(nproc, ntemps), testing=TESTING)
