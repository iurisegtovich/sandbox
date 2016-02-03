import os,sys
from scipy import loadtxt
import matplotlib.pyplot as plt

scores_BACE = loadtxt("GMRQ_scores_BACE_MacroMSMs.txt")
scores_mirco = loadtxt("GMRQ_vs_n_microstates/GMRQ_scores_for_mircoMSMs.txt")
print scores_BACE[:,1]
print scores_mirco
plt.figure()
plt.hold(True)
plt.plot(scores_BACE[:,0],scores_BACE[:,1],label="BACE_macrostates")
plt.plot(scores_mirco[:,0],scores_mirco[:,1],label="Kmeans_mircostates")
plt.ylabel("GMRQ score")
plt.xlabel("Number of states")
plt.title("GMRQ scores of macrostates(BACE) vs. microstates(kmeans)")
plt.legend(loc='best')
plt.savefig("GMRQ_BACE_macro_vs_kmeans_micro.pdf")
plt.show()

