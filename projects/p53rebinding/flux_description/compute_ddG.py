import numpy as np

h_new = 64

ddG = []
binding_affinity = np.loadtxt("binding_affinity_corrected_p53(mM)_vs_mdm2.txt")
binding_affinity_new_helicity = np.loadtxt("binding_affinity_corrected_p53(mM)_vs_mdm2_helicity%d.txt"%h_new)

for i in range(len(binding_affinity)):
    dP =binding_affinity_new_helicity[i]/binding_affinity[i]
    ddG.append(-8.314*300*np.log(dP)*0.001*0.239)
print ddG
output_fn = "ddG_%d_h_corrected_rates.txt"%(h_new)
np.savetxt(output_fn,ddG)
print "Saved: %s"%output_fn
