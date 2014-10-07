from mdtraj import io
import numpy as np

Nv = io.loadh('Nv.h5','arr_0')
tau = 100
Dynamic_quantity = []

for i in range(Nv.shape[0]):
    j = 0
    while j+tau < Nv.shape[1] and Nv[i][j+tau] != -1 :
        Dynamic_quantity.append(np.abs(Nv[i][j]-Nv[i][j+tau]))
        j = j+1

print np.mean(Dynamic_quantity)
print np.std(Dynamic_quantity)
