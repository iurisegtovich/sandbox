import os,sys

for pid in range(6383,6391):
    for stateid in range(40):
        cmd = "pymol -c save_images_pymol.py -- %d %d"%(pid,stateid)
        os.system(cmd)
