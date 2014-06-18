#!/usr/bin/python
#TODO: Save the cartoon structure
import __main__
__main__.pymol_argv = [ 'pymol', '-qc']
import pymol
pymol.finish_launching()
fn = open('./generators/GenId-StateIdpair.txt','r')
pdbstatelist = fn.readlines()
fn.close()
for pdb_state in pdbstatelist:
    pdb_file,state_name = pdb_state.split()
    pdb_name = pdb_file.split('.')[0]
    pymol.cmd.load('./generators/%s'%pdb_file, pdb_name)
    pymol.cmd.disable("all")
    pymol.cmd.enable(pdb_name)
    pymol.cmd.png("./generators/stateimages/%s.png"%state_name)
