import os,sys
import numpy as np
import matplotlib.pyplot as plt

def case1():
    verbose=True
    mapid = 40
    pids = range(6383,6391)
    print pids
    mut = ['EEE','EER','ERE','ERR','REE','RER','RRE','RRR']
    for index,p_id in enumerate(pids):
        data_dir = "Data-{}-macro{}".format(p_id,mapid)
        timescales_fn = os.path.join(data_dir,'ImpliedTimescales-10fold.npy')
        msmts = np.load(timescales_fn)
        msmts_mean = msmts.mean(0)
        msmts_std = msmts.std(0)

        plt.figure()
        for rank in range(msmts.shape[0]):
            plt.errorbar([5]*len(msmts_mean),msmts_mean,yerr=msmts_std,fmt="ro")

        plt.yscale('log')
        plt.xlabel('Lag time (0.1 ns)')
        plt.ylabel('Implied time scales (0.1 ns)')
        plt.title("Fs-{}".format(mut[index]))
        fig_fn = 'ImpliedTimescales-{}-Macro{}_with_errorbar.eps'.format(p_id,mapid)
        plt.savefig(fig_fn)
        if verbose:
            print "Wrote: %s"%fig_fn
    #    plt.show()
        plt.close()

def case2_micro_combined():
    verbose=True
    timescales_fn = os.path.join('Data-combined',"ImpliedTimescales-10fold.npy")
    msmts = np.load(timescales_fn)
    msmts_mean = msmts.mean(0)
    msmts_std = msmts.std(0)

    plt.figure()
    for rank in range(msmts.shape[0]):
        plt.errorbar([5]*len(msmts_mean),msmts_mean,yerr=msmts_std,fmt="ro")

    plt.yscale('log')
    plt.xlabel('Lag time (0.1 ns)')
    plt.ylabel('Implied time scales (0.1 ns)')
    plt.title("Fs-combined")
    fig_fn = 'ImpliedTimescales-micro-combined_with_errorbar.eps'
    plt.savefig(fig_fn)
    if verbose:
        print "Wrote: %s"%fig_fn
#    plt.show()
    plt.close()

def case3_macro_combined():
    verbose=True
    timescales_fn = os.path.join('Data-combined-macro40',"ImpliedTimescales-10fold.npy")
    msmts = np.load(timescales_fn)
    msmts_mean = msmts.mean(0)
    msmts_std = msmts.std(0)

    plt.figure()
    for rank in range(msmts.shape[0]):
        plt.errorbar([5]*len(msmts_mean),msmts_mean,yerr=msmts_std,fmt="ro")

    plt.yscale('log')
    plt.xlabel('Lag time (0.1 ns)')
    plt.ylabel('Implied time scales (0.1 ns)')
    plt.title("Fs-combined")
    fig_fn = 'ImpliedTimescales-combined-Macro40_with_errorbar.eps'
    plt.savefig(fig_fn)
    if verbose:
        print "Wrote: %s"%fig_fn
#    plt.show()
    plt.close()

def case4_print_msmts():
    verbose=True
    mapid = 40
    pids = range(6383,6391)
    print pids
    mut = ['EEE','EER','ERE','ERR','REE','RER','RRE','RRR']
    for index,p_id in enumerate(pids):
        data_dir = "Data-{}-macro{}".format(p_id,mapid)
        timescales_fn = os.path.join(data_dir,'ImpliedTimescales-10fold.npy')
        msmts = np.load(timescales_fn)
        msmts_mean = msmts.mean(0)
        msmts_std = msmts.std(0)
        print mut[index],msmts_mean[0]/10,msmts_std[0]/10
    #microstates combined
    timescales_fn = os.path.join('Data-combined',"ImpliedTimescales-10fold.npy")
    msmts = np.load(timescales_fn)
    msmts_mean = msmts.mean(0)
    msmts_std = msmts.std(0)
    print "micro combined",msmts_mean[0]/10,msmts_std[0]/10
    #macrostates combined
    timescales_fn = os.path.join('Data-combined-macro40',"ImpliedTimescales-10fold.npy")
    msmts = np.load(timescales_fn)
    msmts_mean = msmts.mean(0)
    msmts_std = msmts.std(0)
    print "macro combined",msmts_mean[0]/10,msmts_std[0]/10

def case5_plot_all_msmts():
    verbose=True
    mapid = 40
    pids = range(6383,6391)
    mut = ['EEE','EER','ERE','ERR','REE','RER','RRE','RRR']
    all_msmts =[]
    all_stds = []
    for index,p_id in enumerate(pids):
        data_dir = "Data-{}-macro{}".format(p_id,mapid)
        timescales_fn = os.path.join(data_dir,'ImpliedTimescales-10fold.npy')
        msmts = np.load(timescales_fn)
        msmts_mean = msmts.mean(0)
        msmts_std = msmts.std(0)
        all_msmts.append(msmts_mean/10)
        all_stds.append(msmts_std/10)
    #microstates combined
    timescales_fn = os.path.join('Data-combined',"ImpliedTimescales-10fold.npy")
    msmts = np.load(timescales_fn)
    msmts_mean = msmts.mean(0)
    msmts_std = msmts.std(0)
    all_msmts.append(msmts_mean/10)
    all_stds.append(msmts_std/10)
    #macrostates combined
    timescales_fn = os.path.join('Data-combined-macro40',"ImpliedTimescales-10fold.npy")
    msmts = np.load(timescales_fn)
    msmts_mean = msmts.mean(0)
    msmts_std = msmts.std(0)
    all_msmts.append(msmts_mean/10)
    all_stds.append(msmts_std/10)

    label = mut+ ['Micro','Macro']
    plt.figure()
    plt.hold(True)
    for rank in range(len(all_msmts[0])):
        plt.errorbar(range(len(all_msmts)),np.array(all_msmts)[:,rank],yerr=np.array(all_stds)[:,rank],fmt='o--')
    plt.xticks(range(len(all_msmts)),label)
    plt.ylabel("Implied timescales(ns)")
    plt.savefig("ImpliedTimescales-macro40-all.pdf")


if __name__ == "__main__" :
    #case2_micro_combined()
    #case3_macro_combined()
    #case4_print_msmts()
    case5_plot_all_msmts()
