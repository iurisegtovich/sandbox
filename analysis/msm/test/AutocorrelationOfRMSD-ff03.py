import os,sys
from scipy import loadtxt
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

from scipy.optimize import curve_fit

from msmbuilder import Serializer


RUN = True # If Ture, run the cmd, False, print cmd
DEBUG = False
TESTING = DEBUG

def run_cmd(cmd,testing=False):

    """ If true,RUN, run the cmd, if False print the cmd"""

    if testing:
        os.system(cmd)
    else:
        print '>>',cmd
	
def exp_decay(x, tau0):
    return np.exp(-(x/tau0)) 

	    
def get_RMSD():
    
    rmsds = Serializer.LoadFromHDF('RMSDCluster4.2/Data/RMSD-pdb-gen0.h5')
   
    return rmsds['Data']

def get_projectinfo():

    projectinfo = Serializer.LoadFromHDF('RMSDCluster4.2/ProjectInfo.h5')
    
    return projectinfo



def cal_autocorrelation(traj,ntau):
    
    """calculate the autocorrelation """
    
    trajtrunc = traj[:ntau]
    g = []
    g_trials = []
    if (0): # The SLOW way
	for tau in range(ntau):
	    g.append( np.array([traj[i]*traj[i+tau] for i in range(0, ntau - tau)]).mean() )
    else: # The FAST way
        result = np.correlate(trajtrunc, trajtrunc, mode="full", old_behavior=True)
        print 'len(result)', len(result)
        print result[ntau-1::-1]
        print result[ntau-1:]
	g = result[ntau-1:]
	n = np.arange(ntau,0,-1)
    g_trials = g/n
    return g_trials 
    
    
    
    
	    
    
    

#-----------------------------------------
#------------------MAIN-------------------
#-----------------------------------------
#try:

RMSDs = get_RMSD()
ProjectInfo = get_projectinfo()
Trajlengths = ProjectInfo['TrajLengths']
print 'RMSDs:',RMSDs
print 'Trajlengths:',Trajlengths

g_trials = []
g_mean = []
g_uncertainty = []
time = []
#for i in range(len(alltraj_nhelix_time)):
#    time.append(len(alltraj_nhelix_time[i]))
#minitime = min(time)
#maxtime = max(time)

maxtau = 6000

for i in range(len(RMSDs)):
    g = cal_autocorrelation(RMSDs[i],maxtau) #min(Trajlengths)
    g_trials.append(g)
print g_trials


g_trials_ma = np.ma.empty((len(RMSDs),max(Trajlengths)))
g_trials_ma.mask = True
for i in range(len(g_trials)):
    g_trials_ma[i,:len(g_trials[i])] = g_trials[i]
   
print 'g_trials_ma',g_trials_ma.shape

g_mean = g_trials_ma.mean(axis=0)  
g_uncertainty = g_trials_ma.std(axis=0)[0:maxtau]
#g_uncertainty = g_trials_ma.std(axis=0)/np.sqrt(len(RMSDs))  # uncertainty is calculated over the trials  
#g_uncertainty = g_trials_ma.std(axis=0)[0:maxtau]/np.sqrt(maxtau - np.arange(maxtau)- 1)
#g_uncertainty = g_trials_ma.std(axis=0)[0:maxtau]*(np.arange(maxtau)+1)

#trajcat = []
#for i in range(len(alltraj_nhelix_time)):
#    trajcat = trajcat + alltraj_nhelix_time[i]
    #trajcat = trajcat + alltraj_nhelix_time[i][0:minitime]


# scale the autocorrelation to [0,1]
g = g_mean - (g_trials_ma.mean())**2 # subtract the uncorrelated basline
print 'g',len(g)
g_uncertainty = g_uncertainty/g[0] # then, normalize
print 'g_uncertainty', g_uncertainty
g = g/g[0]
"""g_uncertainty[0:50] *= 0.5
g_uncertainty[50:100] *= 1
g_uncertainty[100:200] *= 2
g_uncertainty[200:500] *= 5
g_uncertainty[500:1000] *= 10
g_uncertainty[1000:maxtau] *= 50
g_uncertainty[maxtau:3000] *= 100
g_uncertainty[3000:4000] *= 200
g_uncertainty[4000:] *= 300
"""

# fit the autocorrelation to a single-exponential curve
v0 = [1000.0]  # Initial guess [tau] for exp(-(x/tau0))
time = np.arange(maxtau)
popt, pcov = curve_fit(exp_decay, time, g[0:maxtau], sigma=g_uncertainty, p0=v0, maxfev=1000000)  # ignore last bin, which has 0 counts
yFit_data = exp_decay(time, popt[0])
print 'best-fit tau_0 = ', popt[0], '+/-', pcov[0][0]


plt.figure()
plt.semilogx(range(maxtau), g[0:maxtau], 'b-')
plt.semilogx(range(maxtau), g[0:maxtau]+g_uncertainty, 'k--',linewidth=0.5)
plt.semilogx(range(maxtau), g[0:maxtau]-g_uncertainty, 'k--',linewidth=0.5)
# plot the fit curve
plt.hold(True)
plt.semilogx(time, yFit_data, 'r-', linewidth=1.0)

plt.xlabel('$\\tau$ (number of steps)')
plt.ylabel('G($\\tau$)')
plt.title('Autocorrelation')

plt.axis([1, maxtau, -0.5, 1.5])
plt.text( np.array(Trajlengths).mean()/100, 1.0, '$\\tau_0$ = %5.2f +/- %5.2f'%(popt[0], pcov[0][0]) )

if (1):
    plt.show()
elif (0):
    #outpdf = 'my-awesome.pdf'
    #print 'Saving PDF of the plot to', outpdf, '...'

    outpng = 'autocorrelation-RMSD-gen0.png'
    print 'Saving PNG of the plot to', outpng, '...'
    plt.savefig(outpng)
    




#except:
    #print "Commands didn't work"
