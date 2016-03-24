from matplotlib import pyplot as plt
from scipy.io import  mmread

from estimator import *

try:
    # msmb 3.0
    from msmbuilder.msm.core import *

except:
    # msmb 2.7
    from msmbuilder import msm_analysis as msm


def get_equilibrium_distribution(T, lagtime=10.0):
    """Given an MSM transition matrix with elements T_ij denoting transition rates from i to j,
    return the normalized stationary eigenvector, and implied timescales.

    INPUT
    T - A transition matrix with elements T_ij denoting transition rates from i to j

    OPTIONS
    lagtime   -  for calculating implied timescales, in units of nanoseconds

    RETURNS
    equil_populations - the normalized stationary eigenvector
    imp_times         - the implied timescales in nanoseconds
    """
    
    nstates = T.shape[0]
    
    # Calculate stationary probability for
    neigs = min(10,nstates)
    #try:
        # msmb 3.0
    evals, right, left = _solve_msm_eigensystem(T,neigs)

    #except:
        # msmb 2.7
    #    evals, right = msm.get_eigenvectors(T, neigs, epsilon=0.001, dense_cutoff=50, right=False, tol=1e-30) #, normalized=True)
    #    evals, left = msm.get_eigenvectors(T, neigs, epsilon=0.001, dense_cutoff=50, right=True, tol=1e-30) #, normalized=True)

    # get the equilibrium populations as the stationary eigenvector
    pi = right[:,0]
    equil_populations = pi/pi.sum()

    # get the implied timescales
    imp_times =  -1.0*lagtime/np.log(evals[1:])
    
    # return normalized populations
    return equil_populations, imp_times


msmts_mut_over_wt = []
msmts_pred_over_wt = []
    

for wt_pid in range(6383,6391):
    tProb_wt_fn = 'Fs-%d-macro40-tProb.mtx'%wt_pid
    T_star = mmread(tProb_wt_fn).todense()
    pi_star, imp_times_star = get_equilibrium_distribution(T_star)

    T_star_flat = np.array(T_star).flatten()

    for mut_pid in range(6383,6391):
        if mut_pid == wt_pid:
            continue

        # Read in the perturbed 150-macrostate transition matrix
        tProb_mut_fn = 'Fs-%d-macro40-tProb.mtx'%mut_pid
        T_perturbed = mmread(tProb_mut_fn).todense()

        pi, imp_times_perturbed = get_equilibrium_distribution(T_perturbed)
    
        # estimate the new rates
        T = estimate_perturbed_rates(T_star, pi_star, pi)

        # get the new implied timescales
        eq, imp_times  = get_equilibrium_distribution(T)

        msmts_mut_over_wt.append(imp_times_perturbed[0]/imp_times_star[0])
        msmts_pred_over_wt.append(imp_times[0]/imp_times_star[0])





    #############################################
plt.figure()
plt.plot(msmts_mut_over_wt, msmts_pred_over_wt, '.')
plt.plot([0.1,10],[0.1,10],'k-')


#R = np.corrcoef(np.log(T_star_flat[Ind]), np.log(T_perturbed_flat[Ind]))[0,1]
#print '$R^2(\\ln p_{ij})$ = %2.3f'%R**2.0
#plt.text(4e-7,1e-1,'$R^2(\\ln p_{ij})$ = %2.3f'%R**2.0)

# calculate rmsd deviation of log p_ij values
#devs = (np.log(T_star_flat[Ind]) -  np.log(T_perturbed_flat[Ind]))
#rmsd = np.sqrt(np.dot(devs,devs)/len(devs))
#print 'rmsd$(\\ln p_{ij})$ = %2.3f'%rmsd
#plt.text(4e-7,1e-2,'rmsd$(\\ln p_{ij})$ = %2.3f'%rmsd)

plt.xscale('log')
plt.yscale('log')
plt.xlabel(r'$log_{10}(\tau_{mut}/\tau_{wt})$')
plt.ylabel(r'$log_{10}(\tau_{pred}/\tau_{wt})$')
plt.savefig('timescale_predictions.pdf')

