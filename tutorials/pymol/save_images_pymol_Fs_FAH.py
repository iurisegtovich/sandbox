"""
#This is a pymol scripts
#To run it, in pymol: @filename

alignto State1-0
show cartoon
hide lines, !State1-0 
center
set cartoon_transparency, 0.5, !State1-0
#color grey70, all
util.chainbow("all")
util.cnc all
"""

#python script
#Run this script with no gui, pymol -c thisfile -- p_id stateid
import os,sys
from pymol import cmd

p_id = int(sys.argv[1])
stateid = int(sys.argv[2])
ref_state_id = 0

outdir=os.path.join("Structure-images","p%d"%p_id)

if not os.path.exists(outdir):
    os.makedirs(outdir)

for i in range(5):
    cmd.load("./PDBs-%d/State%d-%d.pdb"%(p_id,stateid,i))

for i in range(5):
    cmd.align("State%d-%d"%(stateid,i),"State%d-%d"%(stateid,ref_state_id))

#alignto does the same thing as the above loop, but it sometimes could cause problem
#cmd.alignto("State%d-0"%stateid)
cmd.show("cartoon")
cmd.hide("lines","!State%d-%d"%(stateid,ref_state_id))
cmd.center("all")
cmd.set("cartoon_transparency",0.7,"all")
cmd.set("cartoon_transparency",0,"State%d-%d"%(stateid,ref_state_id))
cmd.util.chainbow("all")
cmd.util.cnc("all")
cmd.orient()
cmd.util.performance(0)
cmd.set("ray_opaque_background","off")
cmd.ray()
png_fn = os.path.join(outdir,"image-p%d-state-%d.png"%(p_id,stateid))
high_res_state = [0,1,2,4,11,12,13,20,28,30]
if stateid in high_res_state:
    png_fn = os.path.join(outdir,"image-p%d-state-%d-highres.png"%(p_id,stateid))
    cmd.ray(2400,1800)
    cmd.png(png_fn,dpi=300)
else:
    cmd.png(png_fn)
cmd.quit()

