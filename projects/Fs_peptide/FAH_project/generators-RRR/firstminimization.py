import os, sys, glob, string, commands


def runcmd(cmd,debug):
    if not debug:
        print cmd
        os.system(cmd)
    else:
        print cmd

debug=False

#for gen in range(ngens):
for gen in range(100):
    gro_minimized_hbonds = "Gen%d_solvent_ions_minimized1.gro"%gen
    orientation = 0.0
    fail = False
    while (not os.path.exists(gro_minimized_hbonds)) or fail:
       
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

#structure 
#  number 1 
#  inside cube 0. 0. 0. 45. 
#end structure 

structure %s
  number 1
  fixed 22.5 22.5 22.5 %1f %1f %1f 
  centerofmass
end structure

structure water.pdb
  number 2900 
  inside cube 0. 0. 0. 45.
end structure"""%(gen,starting_pdbfile,orientation,orientation,orientation))
        fout.close()
        orientation += 0.1 
        # Run packmol to get a tmp.pdb with hairpin + exactly 1000 water molecules in a 32A cubic box
        runcmd('packmol < packmol.inp',debug)
    
        # Convert this to a grofile with a cubic box (needed for SMP!) 
        #     make box boundaries 0.5 nm from protein
        #     (the protein will be centered in the box automatically)
        runcmd('editconf -f tmp%d.pdb -o Gen%d_solvated_starting.gro -bt cubic -box 4.5 4.5 4.5'%(gen,gen),debug)
    
        # make a new topology file
        runcmd('echo "6\n1\n" | pdb2gmx -f Gen%d_solvated_starting.gro -o Gen%d_solvated.gro -p Gen%d_solvated.top -ignh'%(gen, gen, gen),debug)
        # make an index file
        runcmd('echo "q\n" | make_ndx -f Gen%d_solvated.gro -o index.ndx'%gen,debug)
    
        # Add counterions at ~100 mM to neutralize the system
        runcmd('grompp -f mdp/minimize1.mdp -c Gen%d_solvated.gro -p Gen%d_solvated.top -n index.ndx -o minimize.tpr'%(gen,gen),debug)
        runcmd('echo "13\n" | genion -s minimize.tpr -n index.ndx -p Gen%d_solvated.top -o Gen%d_solvent_ions.gro -pname NA -nname CL -neutral -conc 0.1'%(gen,gen),debug)
        runcmd('echo "q\n" | make_ndx -f Gen%d_solvent_ions.gro -o index.ndx'%gen,debug)
    
        # Minimize 1
        runcmd('grompp -f mdp/minimize1.mdp -c Gen%d_solvent_ions.gro -p Gen%d_solvated.top -n index.ndx -o minimize2.tpr'%(gen,gen),debug)
        runcmd('mdrun -v -nt 8 -s minimize2.tpr -c Gen%d_solvent_ions_minimized1.gro -g Gen%d_minimize1.log'%(gen,gen),debug)

        # Check the maximum force in log file.
        try:
            force = commands.getoutput('cat Gen%d_minimize1.log | grep "Maximum force"'%gen).split()[3]
            print "Maximum force:",force
            try:
                if float(force)>=100.0:
                    fail = True
                else:
                    fail = False
            except:
                fail = True
        except:
            fail = True
        runcmd('rm \#* *.tpr',debug)
