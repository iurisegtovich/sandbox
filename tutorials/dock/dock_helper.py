import os,sys
from sklearn.externals.joblib import Parallel, delayed

grid_in = ""
rec_path = ""
sphere = ""
dock_in = ""
rec_box = ""


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
            os.system("sed -i -e '/s/rec_struct.mol2/%s/g' %s"%(rec_fn,grid_in)) # modify grid.in file
            #os.system("sed -i -e '/s/state0-0-UNAG/%s/g' %s"%())

def gen_paths(n_jobs,total_n_paths,base_path):
    n_paths_per_job = total_n_paths/n_jobs
    n_paths_left = n_paths_per_job
    c = 0
    all_args = []
    args_per_job = []
    for i in range(100):
        for j in range(5):
            new_working_dir = os.path.join(base_path, "state%d-%d" % (i, j))
            if c < n_paths_per_job:
                args_per_job.append(new_working_dir)
            elif c == n_paths_per_job or c == n_paths_left:
                all_args.append(args_per_job)
                args_per_job = []
                n_paths_left -= c
                c = 0
            c += 1

    return all_args



def run(paths):
    for path in paths:
        os.chdir(paths)
        os.system("grid -i grid.in")
        os.system("dock6 -i anchor_grow_dock.in")


n_jobs = 4
total_n_paths = 500
base_path = ""

setup_working_dirs(base_path)

indices_args = gen_paths(n_jobs,total_n_paths,base_path)
print len(indices_args)

Parallel(n_jobs=n_jobs,verbose=True)(delayed(run)(indices) for indices in indices_args)

