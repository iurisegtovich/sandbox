import os,sys
import seaborn as sns
import numpy as np
from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
from matplotlib import cm

Z = np.array([4,5,6,7])
X = np.array([1,1,2,2])
Y = np.array([1,2,1,2])

#sns.jointplot(x,y,xlim=(x.min(),x.max()),ylim=(y.min(),y.max()),kind='kde')
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(X, Y, Z )

plt.show()
