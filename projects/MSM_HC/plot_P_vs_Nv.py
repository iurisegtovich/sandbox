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

plt.figure()
plt.plot(range(22),p/p.sum(),'o')
plt.plot(range(22),p/p.sum(),'k--')
plt.xlabel('Number of helical residues')
plt.ylabel('Population')
plt.title('Population vs Number of helical residues')
plt.savefig('P_vs_Nv.png',dpi=700)
plt.show()
    
    
