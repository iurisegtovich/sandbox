import os,sys,commands

def run_cmd(cmd)
    print '>>', cmd
    os.system(cmd)


def write_qsub(walltime, run, round, queue, ppn, nthreads,previousJOBID,previous_tpr,extend_by,next_tpr):
    if previousJOBID is not None:
        dependency = '-W depend=afterok:%s'%previousJOBID
    else:
        dependency = ''
    qsubfile = 'qsub%d-ext%d.sh'%(run,round)
    fout = open(qsubfile,'w')
    fout.write("""#!/bin/sh
#PBS -l walltime=%s
#PBS -N RUN%d-ext%d
#PBS -q %s
#PBS -l nodes=1:ppn=%d
#PBS -o RUN%d-ext%d
#PBS %s

cd $PBS_O_WORKDIR

module load gromacs
module load openmpi
"""%(walltime, run, round, queue, ppn, run, round, dependency))
    fout.write("tpbconv -s %s -extend %d -o %s\n"%(previous_tpr,extend_by,next_tpr))
    fout.write("mdrun -v -nt %d -s %s -cpi state%d.cpt -cpo state%d.cpt -noappend\n"%(nthreads,next_tpr,round-1, round))
    #fout.write("time mpirun mdrun_mpi -s %s -cpi state%d.cpt -cpo state%d.cpt -noappend\n"%(next_tpr, round-1, round ))
    fout.close()


extend_by = 20000 # in ps
walltime = '12:00:00'
ppn = 12  # for normal queue
queue = 'normal'
nthreads = 1

run = int(os.path.basename(os.path.abspath(os.curdir)).replace('RUN',''))

if len(sys.argv) < 3:
    print 'Usage: python extend_jobchain.py startround lastround'
    sys.exit(1)
startround = int(sys.argv[1])
lastround = int(sys.argv[2])

if (1):
    previousJOBID = None  # note this will be string when filled
    rounds = range(startround, lastround+1)
    for round in rounds:
        if round == 2:
            os.system('cp production.tpr next1.tpr')
            os.system('cp state_prev.cpt state1.cpt')  # ONLY for a quick fix!!!!
            previous_tpr = 'next%d.tpr'%(round-1)
            next_tpr = 'next%d.tpr'%round
            write_qsub(walltime, run, round, queue, ppn, nthreads,previousJOBID,previous_tpr,extend_by,next_tpr)
            #qsubout = commands.getoutput('qsub %s'%qsubfile)
            #previousJOBID = qsubout.split('.')[0]
            #print 'previousJOBID', previousJOBID




