import os,sys
import numpy as np
import matplotlib.pyplot as plt

verbose=True
mapid = 40
pids = ['combined']+range(6383,6391)
print pids
mut = ['Combined','EEE','EER','ERE','ERR','REE','RER','RRE','RRR']
for index,p_id in enumerate(pids):
    data_dir = "Data-{}-macro{}".format(p_id,mapid)
    lagtime_fn = os.path.join(data_dir,'lagtimes.txt')
    timescales_fn = os.path.join(data_dir,'ImpliedTimescales.npy')
    lagtimes = np.loadtxt(lagtime_fn)
    msmts = np.load(timescales_fn)

    plt.figure()
    for rank in range(msmts.shape[0]):
        plt.plot(lagtimes,[msmts[i][rank] for i in range(len(lagtimes))])

    plt.yscale('log')
    plt.xlabel('Lag time (0.1 ns)')
    plt.ylabel('Implied time scales (0.1 ns)')
    plt.title("Fs-{}".format(mut[index]))
    fig_fn = 'ImpliedTimescales-{}-Macro{}.pdf'.format(p_id,mapid)
    plt.savefig(fig_fn)
    if verbose:
        print "Wrote: %s"%fig_fn
    #plt.show()
    plt.close()


