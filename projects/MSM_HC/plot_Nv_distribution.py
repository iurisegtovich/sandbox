from mdtraj import io
import matplotlib.pyplot as plt
from collections import Counter
import numpy as np

a = io.loadh('results/Nv.h5','arr_0')
a = a.reshape(1,-1)
c = Counter(a[0])
p = np.zeros(22)
for i in range(22):
    p[i] = c[i]

p_MCsampling = np.loadtxt('Nv_distribution_MCSampling.dat')

plt.figure()
plt.plot(range(22),p/p.sum(),'o--')
plt.plot(range(22),p_MCsampling,'o--')
plt.xlabel('Number of helical residues')
plt.ylabel('Population')
plt.title('Population vs Number of helical residues')
plt.legend(['MD simulation','MCM sampling'],prop={'size':10})
plt.savefig('Nv_distribution.png',dpi=700)
#plt.show()
    
    
