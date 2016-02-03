import os,sys
import numpy as np
import matplotlib.pyplot as plt

mapid = 40
seq_id = ['EEE','EER','ERE','ERR','REE','RER','RRE','RRR']
for seqid,pid in enumerate(range(6383,6391)):
    data_dir = "Data-{}-macro{}".format(pid,mapid)
    free_energies = np.loadtxt(os.path.join(data_dir,"Free_energies.dat"))
    Nhs_state = np.loadtxt(os.path.join(data_dir,"Nhs_state.dat"))
    top_paths = np.load(os.path.join(data_dir,"top_paths_substract.npy"))
    top_fluxes = np.load(os.path.join(data_dir,"top_fluxes_substract.npy"))
    plt.figure()
    for i in range(mapid):
        plt.plot(Nhs_state[i],free_energies[i],'ro')
        plt.text(Nhs_state[i],free_energies[i],'%d'%i)
    for path in top_paths:
        plt.plot(Nhs_state[path],free_energies[path])
    plt.plot(Nhs_state[top_paths[0]],free_energies[top_paths[0]],lw=3,c='k')
    print pid
    print top_paths
    print top_fluxes

    fig_fn = "Free_energies_vs_Nhs_state_p{}.pdf".format(pid)
    plt.ylabel("Free energy (kcal/mol)")
    plt.xlabel("Average number of helical residues")
    plt.title("Fs-{}".format(seq_id[seqid]))
#    plt.savefig(fig_fn)
#    print "Saved: %s"%fig_fn
#    plt.show()

        
    
