import os,sys
import numpy as np

#k = 1.3806488*10**-23 #J/K
R = 1.9872041*10**-3 #kcal/(K*mol)
T = 300 #K
for pid in range(6383,6391):
    free_energies_5fold = []
    data_dir = "Data-{pid}-macro40".format(pid=pid)
    Pop_fn = os.path.join(data_dir,"Populations-10fold.npy")
    population_5fold = np.load(Pop_fn)
    for population in population_5fold:
        free_energies = []
        for i in range(len(population)):
            free_energies.append(-R*T*np.log(population[i]))
        free_energies = free_energies - min(free_energies)
        free_energies_5fold.append(free_energies)
    output_fn = os.path.join(data_dir,"Free_energies_10fold.npy")
    np.save(output_fn,free_energies_5fold)
