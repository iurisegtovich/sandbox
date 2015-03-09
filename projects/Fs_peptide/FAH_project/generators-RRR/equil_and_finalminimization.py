import os, sys, glob, string, commands


def runcmd(cmd,debug):
    if not debug:
        print cmd
        os.system(cmd)
    else:
        print cmd

debug=False

#for gen in range(ngens):
for gen in range(6,10):

    Gro_minimized = "Gen%d_solvent_ions_minimized2.gro"%gen
    if os.path.exists(Gro_minimized):
        # Equilibrate with pressure coupling on (Berendsen)
        runcmd('grompp -f mdp/equil.mdp -c %s -p Gen%d_solvated.top -n index.ndx -o equil.tpr'%(Gro_minimized,gen),debug)
        runcmd('mdrun -v -nt 8 -s equil.tpr -c Gen%d_solvent_ions_equilibrated.gro -e ener-equil%d.edr -g equil%d.log -o traj-equil%d.trr -x traj-equil%d.xtc' %(gen,gen,gen,gen,gen),debug)
    
        # Minimize one last time (to be safe, as FAH will assign random velocities)
        runcmd('grompp -f mdp/minimize_final.mdp -c Gen%d_solvent_ions_equilibrated.gro -p Gen%d_solvated.top -n index.ndx -o minimize3.tpr'%(gen,gen),debug)
        runcmd('mdrun -v -nt 8 -s minimize3.tpr -c Gen%d_ready4fah.gro -e ener-finalmin%d.edr -g finalmin%d.log -o traj-finalmin%d.trr -x traj-finalmin%d.xtc'%(gen,gen,gen,gen,gen),debug)
    
        # delete backup files and mdp tpr files
        runcmd('rm \#* *.tpr',debug)
    else:
        print "Can't find %s"%Gro_minimized
        system.exit()
