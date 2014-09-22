import numpy as np
import matplotlib.pyplot as plt
import argparse

parser=argparse.ArgumentParser()
parser.add_argument('-p','--Populations',help='Populations file. Default: Populations.dat',default='Populations.dat')
args = parser.parse_args()

p = np.loadtxt(args.Populations)

plt.figure()
plt.plot(range(len(p)),p,'o')
plt.show()