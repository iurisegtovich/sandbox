import os, sys, glob, string

import numpy as np

import matplotlib.pyplot as plt
from scipy import loadtxt

#from msmbuilder import Serializer
from msmbuilder import io
from scipy.io import mmread

from msmbuilder.MSMLib import *

import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

# default settings
# VAV:  For some reason, this needs to be set in the main body of the script, not in a subroutine
plt.rc('figure', figsize=(6.5, 3.3))  # two-col width can be between 4.167 and 7 in
plt.rc('lines', linewidth=1.5, markersize=5)
plt.rc('font', size=16.0)
fontfamily = {'family':'sans-serif','sans-serif':['Arial']}
plt.rc('font', **fontfamily)
plt.rc('xtick', labelsize='small')
plt.rc('ytick', labelsize='small')
plt.rc('legend', fontsize=12.0)

sys.path.append('../')
from Figure import Figure



def plot_landscape(TPTdir, unfolded_state=None, folded_state=None, UseLabels=False, UseNewArrows=True, title=None):

    # Load in equilibrium populations
    print "os.path.join(TPTdir,'*Populations*.dat')", os.path.join(TPTdir,'*Populations*.dat')
    popfile = glob.glob( os.path.join(TPTdir,'*Populations*.dat') )[0]
    print popfile
    p = loadtxt(popfile)

    # Convert populations to free energy
    kT = 0.5959
    F = -kT*np.log(p)
    #F = F - F.min()   # this makes native state F=0
    # F = F - F[unfolded_state]  # this makes the unfolded state F=0
    F = F - F[folded_state]  # this makes the folded state F=0


    # Load in pfold values
    Fcommitfile = glob.glob( os.path.join(TPTdir,'*committors*.dat') )[0]
    pfold = loadtxt(Fcommitfile)

    plt.plot(pfold, F, 'k.') 

    pathfile = glob.glob( os.path.join(TPTdir,'*Paths*.h5') )[0]
    #paths = Serializer.LoadFromHDF(pathfile)
    paths = io.loadh(pathfile)
    toppaths = paths["Paths"]
    print 'toppaths', toppaths
    fluxes = paths['fluxes']
    reminder = """{'Paths': array([[  0.,   1.,   9.,  25.,  67.,  66.,  -1.,  -1.],
       [  0.,   2.,  11.,  15.,  26.,  66.,  -1.,  -1.],
       [  0.,   5.,  25.,  26.,  66.,  -1.,  -1.,  -1.],
       [  0.,   2.,   4.,   8.,  31.,  25.,  26.,  66.],
       [  0.,   2.,   8.,  31.,  25.,  67.,  66.,  -1.],
       [  0.,   2.,  17.,  26.,  66.,  -1.,  -1.,  -1.],
       [  0.,   2.,  15.,  26.,  66.,  -1.,  -1.,  -1.],
       [  0.,   5.,  15.,  26.,  66.,  -1.,  -1.,  -1.],
       [  0.,   2.,  30.,  59.,  66.,  -1.,  -1.,  -1.],
       [  0.,  11.,  43.,  66.,  -1.,  -1.,  -1.,  -1.]]), 'fluxes': array([  3.88076552e-07,   2.99234810e-07,   2.45282787e-07,
         2.05725488e-07,   9.79638218e-08,   9.73667502e-08,
         8.60687117e-08,   7.71866221e-08,   6.33758286e-08,
         6.85657260e-08]), 'Bottlenecks': [], 'SerializerFilename': '/Users/vince/git/voelzlab/surprisal/examples/nina-toymer/TPT-AAFAAA/Paths.h5'}"""

    # compile the total flux for edges seen in the top paths
    totalflux = {}  # make a dict, because not all edges will be in the top paths
    for n in range(toppaths.shape[0]):
        i = 0
        while (i<(toppaths.shape[1]-1)) and (toppaths[n,i+1] >= 0):
             j, k = toppaths[n,i], toppaths[n,i+1]
             if not totalflux.has_key((j,k)):
                 totalflux[(j,k)] = 0.0
             totalflux[(j,k)] += fluxes[i]
             i += 1
    print 'totalflux', totalflux
    c = FluxConverter(np.array(totalflux.values()))   # object for convert flux values in line widths

    
    #surprisals = loadtxt('surprisals.dat')
    #print 'surprisals', surprisals
    #s = FluxConverter(surprisals, max_lw=20.0)
    #s.set_minflux(1e-10)
    #s.set_maxflux(1e-4)

    # plot the nodes
    j = 0
    while j < len(pfold):
        #plt.plot( pfold[j], F[j], 'o', markersize=s.convert(surprisals[j]) )
        if UseLabels:
            plt.text( pfold[j], F[j], str(int(j)), size=10)
        plt.hold(True)
        j += 1

    maxpaths = 10
    #for n in range(1,toppaths.shape[0])+[0]:   # plot top path (n=0) last
    for n in range(1, maxpaths)+[0]:
        print 'plotting path', n, toppaths[n,:]
        # plot the edges
        i = 0
        while (i<(toppaths.shape[1]-1)) and (toppaths[n,i+1] >= 0):
            j,k = toppaths[n,i], toppaths[n,i+1] 
            #plt.plot( [pfold[j], pfold[k]] , [F[j], F[k]], 'k-', linewidth=c.convert(totalflux[(j,k)]) )
            if not UseNewArrows:
                # --- old arrows ---
                arr = plt.Arrow(pfold[j], F[j], (pfold[k]-pfold[j]), (F[k]-F[j]), edgecolor='none', width=0.1*c.convert(totalflux[(j,k)]) )
                ax.add_patch(arr) 
                if n > 0:
                    arr.set_facecolor('k')
                else:
                    arr.set_facecolor('r')
            else:
                # --- new arrows ---  
                if n > 0:
                    arrowcolor = 'k'
                else:
                    arrowcolor = 'r'
                #plt.arrow( x, y, dx, dy, **kwargs )
                plt.arrow( pfold[j], F[j], 0.95*(pfold[k]-pfold[j]), 0.95*(F[k]-F[j]), fc=arrowcolor, ec=arrowcolor,
                       head_width=0.02, head_length=0.03, width=0.010*c.convert(totalflux[(j,k)]) )
            plt.hold(True)
            i += 1

    plt.xlabel('$p_{fold}$')
    plt.ylabel('free energy (kcal/mol)')
    plt.axis( [-0.1, 1.1, -1, 1.5] )
    plt.xticks(np.arange(0,1.1,0.1))


    if title != None:
        plt.title(title)


class FluxConverter(object):

    def __init__(self, fluxes, max_lw=5.0, min_lw=0.3):

        self.minflux = fluxes.min()
        self.maxflux = fluxes.max()
  
        self.max_lw = max_lw
        self.min_lw = min_lw

    def set_minflux(self, minflux):
        self.minflux = minflux

    def set_maxflux(self, maxflux):
        self.maxflux = maxflux

    def convert(self, flux, logScale=True):
        
        if logScale:
            logflux = np.log(flux)
            fraction = max(0.0, min(1.0, (logflux - np.log(self.minflux))/(np.log(self.maxflux) - np.log(self.minflux))))
            return np.exp( np.log(self.min_lw) + fraction* (np.log(self.max_lw) - np.log(self.min_lw)))
        else:
            fraction = max(0.0, min(1.0, (flux - np.log(self.minflux))/(self.maxflux - self.minflux)))
            return min_lw + fraction*(self.max_lw - self.min_lw)
             
        return np.real(evecs[:,0]/evecs[:,0].sum())

def H(row):
    result = 0.
    for j in range(row.shape[0]):
        if row[i] > 0.:
            result += -row[i]*np.log(row[i])
    return result

if __name__ == '__main__':

    usage = """Usage:  python get_populations.py TPTdir1 TPTdir2
        will write outname.Populations.dat"""

    if len(sys.argv) < 3:
        print usage
        sys.exit(1)

    TPTdir1 = sys.argv[1]
    TPTdir2 = sys.argv[2]

    # Create the figure
    f = Figure(plt)
    f.add_panel(0.07, 0.12, 0.42, 0.8)      # [xorig, yorig, width, height]
    f.add_panel(0.57, 0.12, 0.42, 0.8)
   
    ################################
    # Panel A 
    ax = plt.axes(f.axis_handles[0])
    plot_landscape(TPTdir1, unfolded_state=1, folded_state=10, UseLabels=True, title='GB1') 

    ###############################
    # Panel B 
    ax = plt.axes(f.axis_handles[1])
    plot_landscape(TPTdir2, unfolded_state=1, folded_state=10, UseLabels=True, title='trpzip4')


# plt.show()


fig = plt.gcf()
fig.set_size_inches(5.75,2.5)

outfile = 'tpt_wt_vs_mut.eps'
plt.savefig(outfile, format='eps')
print 'Wrote:', outfile


