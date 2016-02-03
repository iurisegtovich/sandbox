import os,sys
import numpy as np

#k = 1.3806488*10**-23 #J/K
R = 1.9872041*10**-3 #kcal/(K*mol)
T = 300 #K
for pid in range(6383,6391):
    free_energies = []
    data_dir = "Data-{pid}-macro40".format(pid=pid)
    Pop_fn = os.path.join(data_dir,"Populations.dat")
    population = np.loadtxt(Pop_fn)
    for i in range(len(population)):
        free_energies.append(-R*T*np.log(population[i]))
    free_energies = free_energies - min(free_energies)
    output_fn = os.path.join(data_dir,"Free_energies.dat")
    np.savetxt(output_fn,free_energies)
