import os,sys
from sklearn.externals.joblib import Parallel, delayed

def setup_working_dirs(base_path):
    #base_path = "~/projects"
    for i in range(100):
        for j in range(5):
            new_working_dir = os.path.join(base_path,"state%d-%d"%(i,j))
            if not os.path.exists(new_working_dir):
                os.makedirs(new_working_dir)
            os.system("cp ./prep/*.in %s"%new_working_dir)
            grid_in = os.path.join(new_working_dir,"grid.in")
            dock_in = os.path.join(new_working_dir,"dock.in")
            rec_fn = "\/Users\/tud51931\/projects\/murA\/MurA-MSM-mol2\/state%d-%d.aligned.mol2"%(i,j)
            cmd = "sed -i -e 's/receptor.mol2/%s/g' %s"%(rec_fn,grid_in)
            print cmd
            os.system(cmd) # modify grid.in file
            #os.system("sed -i -e 's/state0-0-UNAG/%s/g' %s"%())

def gen_paths(n_jobs,total_n_paths,base_path):
    n_paths_per_job = total_n_paths/n_jobs
    n_paths_left = n_paths_per_job
    c = 1
    all_args = []
    args_per_job = []
    for i in range(100):
        for j in range(5):
            new_working_dir = os.path.join(base_path, "state%d-%d" % (i, j))
            if c <= n_paths_per_job:
                args_per_job.append(new_working_dir)
            if c == n_paths_per_job or c == n_paths_left:
                all_args.append(args_per_job)
                args_per_job = []
                n_paths_left -= c
                c = 0
            c += 1

    return all_args



def run_grid(paths):
    for path in paths:
        os.chdir(path)
        os.system("grid -i grid.in")

def run_dock6(paths):
    for path in paths:
        os.chdir(path)
        os.system("dock6 -i anchor_grow_dock.in")


n_jobs = 4
total_n_paths = 500
base_path = "/Users/tud51931/projects/murA/MurA-dock-MSMs"
#rec_path = "/Users/tud51931/projects/murA/MurA-MSM-mol2"

#setup_working_dirs(base_path)

indices_args = gen_paths(n_jobs,total_n_paths,base_path)
#print indices_args
#print len(indices_args)

Parallel(n_jobs=n_jobs,verbose=True)(delayed(run_grid)(indices) for indices in indices_args)
Parallel(n_jobs=n_jobs,verbose=True)(delayed(run_dock6)(indices) for indices in indices_args)

