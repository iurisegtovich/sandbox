import os,sys
import numpy as np
import matplotlib.pyplot as plt

mapid = 40
seq_id = ['EEE','EER','ERE','ERR','REE','RER','RRE','RRR']
plt.figure(figsize=(14,24))
for seqid,pid in enumerate(range(6383,6391)):
    data_dir = "Data-{}-macro{}".format(pid,mapid)
    free_energies = np.loadtxt(os.path.join(data_dir,"Free_energies.dat"))
    Nhs_state = np.loadtxt(os.path.join(data_dir,"Nhs_state.dat"))
    plt.subplot(4,2,seqid+1)
    for i in range(mapid):
        plt.plot(Nhs_state[i],free_energies[i],'ro')
        plt.text(Nhs_state[i],free_energies[i],'%d'%i)
    plt.ylabel("Free energy (kcal/mol)")
    plt.xlabel("Average number of helical residues")
    plt.title("Fs-{}".format(seq_id[seqid]))
#plt.show()
fig_fn = "Free_energies_vs_Nhs_state.pdf".format(pid)
plt.savefig(fig_fn)
print "Saved: %s"%fig_fn

        
    
