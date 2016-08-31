import json
import os,sys
import commands

project_id = 8610
total_clones = 5000
total_runs = 1

project_path = 'PROJ%d'%project_id
converted_path = 'Converted%d'%project_id
logfile_fn = os.path.join(converted_path,'trajectories.jsonl')

if not os.path.exists(converted_path):
    os.makedirs(converted_path)

i = 0
for run in range(total_runs):
    for clone in range(total_clones):
        traj_path = os.path.join(project_path,'RUN%d'%run,'CLONE%d'%clone)
        frame0_path = os.path.join(traj_path,'frame0.xtc')
        if os.path.exists(frame0_path):
            frames_path = os.path.join(traj_path,'frame*.xtc')
            output_traj_fn = 'traj-%08d.xtc'%i 
            output_fn = os.path.join(converted_path,output_traj_fn)
            cmd = 'trjcat -f %s -o %s -overwrite'%(frames_path,output_fn)
            os.system(cmd)
            log_entry = {'filename':output_traj_fn, 'sourcepath':traj_path}
            with open(logfile_fn,'a') as output:
                json.dump(log_entry,output)
                output.write('\n')
            i = i + 1

