import os,sys
import numpy as np
import matplotlib.pyplot as plt

mapid = 40
seq_id = ['EEE','EER','ERE','ERR','REE','RER','RRE','RRR']
plt.figure(figsize=(14,24))
for seqid,pid in enumerate(range(6383,6391)):
    data_dir = "Data-{}-macro{}".format(pid,mapid)
    free_energies = np.load(os.path.join(data_dir,"Free_energies_10fold.npy"))
    free_energies_mean = free_energies.mean(0)
    free_energies_std = free_energies.std(0)
    Nhs_state = np.loadtxt(os.path.join(data_dir,"Nhs_state.dat"))
    top_paths = np.load(os.path.join(data_dir,"top_paths_substract.npy"))
    top_fluxes = np.load(os.path.join(data_dir,"top_fluxes_substract.npy"))
    plt.subplot(4,2,seqid+1)
    for i in range(mapid):
        plt.errorbar(Nhs_state[i],free_energies_mean[i],yerr=free_energies_std[i],fmt='ro')
        plt.text(Nhs_state[i],free_energies_mean[i],'%d'%i)
    for path in top_paths:
        plt.plot(Nhs_state[path],free_energies_mean[path])
    plt.plot(Nhs_state[top_paths[0]],free_energies_mean[top_paths[0]],lw=3,c='k')
    plt.ylabel("Free energy (kcal/mol)")
    plt.xlabel("Average number of helical residues")
    plt.title("Fs-{}".format(seq_id[seqid]))
#plt.show()
fig_fn = "Free_energies_vs_Nhs_state_with_paths_errorbars_together.pdf".format(pid)
plt.savefig(fig_fn)
print "Saved: %s"%fig_fn

        
    
