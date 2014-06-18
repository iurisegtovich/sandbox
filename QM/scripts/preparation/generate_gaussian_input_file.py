import os,sys

#for i in range(99,100):
#    cmd = "echo '\n%%NProcShared=8\n#N B3LYP/6-311+G(2d,p)//HF/6-31G(d)\n\ncineromycineB-%d\n\n0,1\n' | newzmat -prompt -ipdb -ozmat %d.pdb %d.com"%(i,i,i)
    #print cmd
#    os.system(cmd)
input_pdb = 'Nspe5_ACE_NDM.pdb'
output_com = 'Nspe5_ACE_NDM.com'
cmd = "echo '\n%%NProcShared=8\n#N HF/6-31G(d) OPT\n\nNspe5-ACE-NDM\n\n0,1\n' | newzmat -prompt -ipdb -ozmat %s %s"%(input_pdb,output_com)
os.system(cmd)
