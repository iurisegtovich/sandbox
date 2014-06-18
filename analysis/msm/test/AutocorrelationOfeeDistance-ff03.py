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
	
def exp_decay(x, tau0,tau1):
    return np.exp(-(x/tau0))+np.exp(-(x/tau1)) 
    
    
def get_projectinfo():
    
    try:
        projectinfo = Serializer.LoadFromHDF('../ProjectInfo.h5')
    except IOError:
		print "Can't find ProjectInfo.h5!"

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
        result = np.correlate(trajtrunc, trajtrunc, mode="full")
	g = result[ntau-1:]
	n = np.arange(ntau,0,-1)
    g_trials = g/n
    return g_trials 
    
    
    
    
	    
    
    

#-----------------------------------------
#------------------MAIN-------------------
#-----------------------------------------


try:
    Distances = loadtxt('../EndtoEndDistances.dat')
    Distances_masked = np.ma.array(Distances,mask = (Distances == -1))
except IOError:
	print "Can't find EndtoEndDistances.dat"

print 'Distances_masked', Distances_masked,'shape',Distances_masked.shape
ProjectInfo = get_projectinfo()
Trajlengths = ProjectInfo['TrajLengths']

g_trials = []
g_mean = []
g_uncertainty = []
time = []

for i in range(len(Distances_masked)):
    g = cal_autocorrelation(Distances_masked[i],min(Trajlengths)) #min(Trajlengths)
    g_trials.append(g)
print g_trials


g_trials_ma = np.ma.empty((len(Distances_masked),min(Trajlengths)))
g_trials_ma.mask = True
for i in range(len(g_trials)):
    g_trials_ma[i,:len(g_trials[i])] = g_trials[i]
   
print 'g_trials_ma',g_trials_ma

g_mean = g_trials_ma.mean(axis=0)  
#g_uncertainty = g_trials_ma.std(axis=0)[0:maxtau]
g_uncertainty = g_trials_ma.std(axis=0)/np.sqrt(len(Distances_masked))  # uncertainty is calculated over the trials  
print "g_trials_ma.std(axis=0)",g_trials_ma.std(axis=0)
#g_uncertainty = g_trials_ma.std(axis=0)[0:maxtau]/np.sqrt(maxtau - np.arange(maxtau)- 1)
#g_uncertainty = g_trials_ma.std(axis=0)[0:maxtau]*(np.arange(maxtau)+1)

#trajcat = []
#for i in range(len(alltraj_nhelix_time)):
#    trajcat = trajcat + alltraj_nhelix_time[i]
    #trajcat = trajcat + alltraj_nhelix_time[i][0:minitime]


# scale the autocorrelation to [0,1]
g = g_mean - (Distances_masked.mean())**2 # subtract the uncorrelated basline
print 'g_mean',g_mean,"g_trials_ma.mean()",g_trials_ma.mean()
g_uncertainty = g_uncertainty/g[0] # then, normalize
print 'g_uncertainty', g_uncertainty
g = g/g[0]
print 'g',g

# fit the autocorrelation to a single-exponential curve
v0 = [100.0]  # Initial guess [tau] for exp(-(x/tau0))
time = np.arange(len(g))
popt, pcov = curve_fit(exp_decay, time, g, sigma = g_uncertainty, p0=v0, maxfev=10000)  # ignore last bin, which has 0 counts
yFit_data = exp_decay(time, popt[0],popt[1])
print 'best-fit tau_0 = ', popt[0], '+/-', pcov[0][0]


plt.figure()
plt.semilogx(range(len(g)), g[0:len(g)], 'b-')
plt.semilogx(range(len(g)), g[0:len(g)]+g_uncertainty, 'k--',linewidth=0.5)
plt.semilogx(range(len(g)), g[0:len(g)]-g_uncertainty, 'k--',linewidth=0.5)
# plot the fit curve
plt.hold(True)
plt.semilogx(time, yFit_data, 'r-', linewidth=1.0)

plt.xlabel('$\\tau$ (number of steps)')
plt.ylabel('G($\\tau$)')
plt.title('Autocorrelation')

plt.axis([1, len(g), -0.5, 1.5])
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
