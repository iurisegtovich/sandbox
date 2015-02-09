import os, sys, glob, string, commands


def runcmd(cmd,debug):
    if not debug:
        print cmd
        os.system(cmd)
    else:
        print cmd

debug=True

#for gen in range(ngens):
for gen in [2]:

    starting_pdbfile = 'gen%d.pdb'%gen

    # write PDB with a single water in it (for packmol) 
    if not(os.path.exists('water.pdb')):
        fout = open('water.pdb','w')
        fout.write("""ATOM      1  OW  SOL     1      10.060   7.970  17.100  1.00  0.00            
ATOM      2  HW1 SOL     1      10.070   8.550  17.860  1.00  0.00            
ATOM      3  HW2 SOL     1      10.720   8.330  16.510  1.00  0.00            
TER
END
""")
        fout.close()

    # Write the input file for packmol
    # This will get exactly 1500 water molecules in a 36A cubic box 
    fout = open('packmol.inp', 'w')
    fout.write("""tolerance 2.0 
output tmp%d.pdb 
filetype pdb 

structure %s
  number 1 
  inside cube 0. 0. 0. 45. 
end structure 

#structure 
#  number 1
#  fixed 20. 20. 20. 0. 0. 0.
#  centerofmass
#end structure

structure water.pdb
  number 2900 
  inside cube 0. 0. 0. 45.
end structure"""%(gen,starting_pdbfile))
    fout.close()

    # Run packmol to get a tmp.pdb with hairpin + exactly 1000 water molecules in a 32A cubic box
    runcmd('packmol < packmol.inp',debug=True)

    # Convert this to a grofile with a cubic box (needed for SMP!) 
    #     make box boundaries 0.5 nm from protein
    #     (the protein will be centered in the box automatically)
    runcmd('editconf -f tmp%d.pdb -o Gen%d_solvated_starting.gro -bt cubic -box 4.5 4.5 4.5'%(gen,gen),debug)

    # make a new topology file
    runcmd('echo "6\n1\n" | pdb2gmx -f Gen%d_solvated_starting.gro -o Gen%d_solvated.gro -p Gen%d_solvated.top'%(gen, gen, gen),debug)
    # make an index file
    runcmd('echo "q\n" | make_ndx -f Gen%d_solvated.gro -o index.ndx'%gen,debug)

    # Add counterions at ~100 mM to neutralize the system
    runcmd('grompp -f mdp/minimize.mdp -c Gen%d_solvated.gro -p Gen%d_solvated.top -n index.ndx -o minimize.tpr'%(gen,gen),debug)
    runcmd('echo "13\n" | genion -s minimize.tpr -n index.ndx -p Gen%d_solvated.top -o Gen%d_solvent_ions.gro -pname NA -nname CL -neutral -conc 0.1'%(gen,gen),debug)
    runcmd('echo "q\n" | make_ndx -f Gen%d_solvent_ions.gro -o index.ndx'%gen,debug)

    # Minimize
    runcmd('grompp -f mdp/minimize.mdp -c Gen%d_solvent_ions.gro -p Gen%d_solvated.top -n index.ndx -o minimize2.tpr'%(gen,gen),debug)
    runcmd('mdrun -v -nt 8 -s minimize2.tpr -c Gen%d_solvent_ions_minimized.gro'%gen,debug)

    """
  -s      topol.tpr  Input        Run input file: tpr tpb tpa
  -o       traj.trr  Output       Full precision trajectory: trr trj cpt
  -x       traj.xtc  Output, Opt. Compressed trajectory (portable xdr format)
-cpi      state.cpt  Input, Opt.  Checkpoint file
-cpo      state.cpt  Output, Opt. Checkpoint file
  -c    confout.gro  Output       Structure file: gro g96 pdb etc.
  -e       ener.edr  Output       Energy file
  -g         md.log  Output       Log file
    """
    # Equilibrate with pressure coupling on (Berendsen)
    runcmd('grompp -f mdp/equil.mdp -c Gen%d_solvent_ions_minimized.gro -p Gen%d_solvated.top -n index.ndx -o equil.tpr'%(gen,gen),debug)
    runcmd('mdrun -v -nt 8 -s equil.tpr -c Gen%d_solvent_ions_equilibrated.gro -e ener-equil%d.edr -g equil%d.log -o traj-equil%d.trr -x traj-equil%d.xtc' %(gen,gen,gen,gen,gen),debug)

    # Minimize one last time (to be safe, as FAH will assign random velocities)
    runcmd('grompp -f mdp/minimize.mdp -c Gen%d_solvent_ions_equilibrated.gro -p Gen%d_solvated.top -n index.ndx -o minimize3.tpr'%(gen,gen),debug)
    runcmd('mdrun -v -nt 8 -s minimize3.tpr -c Gen%d_ready4fah.gro -e ener-finalmin%d.edr -g finalmin%d.log -o traj-finalmin%d.trr -x traj-finalmin%d.xtc'%(gen,gen,gen,gen,gen),debug)

    # delete backup files and mdp tpr files
    runcmd('rm \#* *.mdp *.tpr',debug)

