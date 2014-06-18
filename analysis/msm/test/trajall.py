import os

RUN = True # If Ture, run the cmd, False, print cmd
DEBUG = False
TESTING = RUN

def runcmd(cmd,testing=False):

    """ If true,RUN, run the cmd, if False print the cmd"""

    if testing:
        os.system(cmd)
    else: 
        print '>>',cmd


for i in range(0,8):
    runcmd('trjcat -f traj%d.xtc traj%d.part000?.xtc -o trajall%d.xtc'%(i,i,i),TESTING)
    runcmd('eneconv -f ener%d.edr ener%d.part000?.edr -o enerall%d.edr'%(i,i,i),TESTING) 
