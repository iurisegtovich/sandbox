import os,sys

#if len(sys.argv) < 2:
#    print "Usage: python convert_whole.py project_id"
#    sys.exit()

project_id = 8610

trajs_path = "Converted%d"%project_id
trajs_path_new = "Converted%d-whole"%project_id
if not os.path.exists(trajs_path_new):
    os.makedirs(trajs_path_new)
json_fn = os.path.join(trajs_path,"trajectories.jsonl")

#tpr = './PROJ%d/RUN0/CLONE0/frame0.tpr'%project_id
tpr = "./p%d/frame0.tpr"%project_id

with open(json_fn,'r') as fn:
    lines = fn.readlines()
    total_trajs = len(lines)

if os.path.exists(tpr):
    for i in range(total_trajs):
        fn_xtc = os.path.join(trajs_path,'traj-%08d.xtc'%i)
        fn_whole_xtc = os.path.join(trajs_path,'traj-%08d.whole.xtc'%i)
        cmd = "echo '1\n' | gmx trjconv -f %s -o %s -pbc whole -s %s"%(fn_xtc,fn_whole_xtc,tpr)
        os.system(cmd)
else:
    print "Cannot find %s"%tpr

