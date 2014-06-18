import os,sys

for i in range(71,100):
    qsub = open('qsub%d.sh'%i,'w')
    content = """#!/bin/sh
#PBS -l walltime=3:00:00
#PBS -N cineromycinB-%d
#PBS -q legacy
#PBS -l nodes=1:ppn=8
#PBS -M tud51931@temple.edu
#PBS
cd $PBS_O_WORKDIR

module load gaussian/g09
g09 -i %d.com
"""%(i,i)
    qsub.write(content)
    qsub.close()

