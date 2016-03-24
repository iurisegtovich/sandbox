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


    
plt.figure(figsize=(8,6))

# Read in the unpertubred 150-macrostate transition matrix for gb1 hairpin 
T_star = mmread('Fs-6383-macro40-tProb.mtx').todense()
pi_star, imp_times_star = get_equilibrium_distribution(T_star)

T_star_flat = np.array(T_star).flatten()

print 'T_star_flat', T_star_flat
print 'T_star_flat.shape', T_star_flat.shape


panel = 1    
identifiers = ['Fs-6390-macro40']
nsystems = len(identifiers)
for identifier in identifiers: #,'tz5','tz6']:

       
    # Read in the perturbed 150-macrostate transition matrix
    T_perturbed = mmread('%s-tProb.mtx'%identifier).todense()
    print 'T_perturbed', T_perturbed

    pi, imp_times_perturbed = get_equilibrium_distribution(T_perturbed)
    
    # estimate the new rates
    T = estimate_perturbed_rates(T_star, pi_star, pi)

    # get the new implied timescales 
    eq, imp_times  = get_equilibrium_distribution(T)
    
    print 'old T', T_star
    print 'final estimate of T:', T


    print 'Rendering scatter plots for', identifier, '...'
    nstates = T.shape[0]
   
    ############################################# 
    plt.subplot(2,2,panel)
    for i in range(nstates):
        for j in range(nstates):
            plt.plot(T_star[i,j], T_perturbed[i,j], '.')

    # calculate an R^2 value for the log p_ij values 
    T_perturbed_flat = np.array(T_perturbed).flatten()
    Ind = (T_star_flat > 0.)*(T_perturbed_flat > 0.)
    R = np.corrcoef(np.log(T_star_flat[Ind]), np.log(T_perturbed_flat[Ind]))[0,1]
    print '$R^2(\\ln p_{ij})$ = %2.3f'%R**2.0
    plt.text(4e-7,1e-1,'$R^2(\\ln p_{ij})$ = %2.3f'%R**2.0)  

    # calculate rmsd deviation of log p_ij values 
    devs = (np.log(T_star_flat[Ind]) -  np.log(T_perturbed_flat[Ind]))
    rmsd = np.sqrt(np.dot(devs,devs)/len(devs))
    print 'rmsd$(\\ln p_{ij})$ = %2.3f'%rmsd
    plt.text(4e-7,1e-2,'rmsd$(\\ln p_{ij})$ = %2.3f'%rmsd)  

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Unperturbed rates (Fs-6383)')
    plt.ylabel('Perturbed rates (%s)'%identifier)
    plt.title(identifier)
    panel += 1


    ##################################
    plt.subplot(2,2,panel)
    for i in range(nstates):
        for j in range(nstates):
            plt.plot(T[i,j], T_perturbed[i,j], '.')

    # calculate an R^2 value for the log p_ij values 
    T_flat = np.array(T).flatten()
    Ind = (T_flat > 0.)*(T_perturbed_flat > 0.)
    R = np.corrcoef(np.log(T_flat[Ind]), np.log(T_perturbed_flat[Ind]))[0,1]
    print '$R^2(\\ln p_{ij})$ = %2.3f'%R**2.0
    plt.text(4e-7,1e-1,'$R^2(\\ln p_{ij})$ = %2.3f'%R**2.0)

    # calculate rmsd deviation of log p_ij values 
    devs = (np.log(T_flat[Ind]) -  np.log(T_perturbed_flat[Ind]))
    rmsd = np.sqrt(np.dot(devs,devs)/len(devs))
    print 'rmsd$(\\ln p_{ij})$ = %2.3f'%rmsd
    plt.text(4e-7,1e-2,'rmsd$(\\ln p_{ij})$ = %2.3f'%rmsd)  

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('RelEntropy estimate (from Fs-6383)')
    plt.ylabel('Perturbed rates (%s)'%identifier)
    plt.title(identifier)
    panel += 1


    ############################################# 
    plt.subplot(2,2,panel)
    # plot implied timescales for gb1 vs perturbed 
    for tau in imp_times_star:
        plt.plot([1,2],[tau,tau], 'k-')
    for tau in imp_times_perturbed:
        plt.plot([3,4],[tau,tau], 'k-')
    plt.xlim(0,5)
    plt.xticks([1.5, 3.5], ['Fs-6383', identifier])
    plt.yscale('log')
    plt.ylabel('implied timescales (ns)')
    panel += 1

    ############################################# 
    plt.subplot(2,2,panel)
    # plot implied timescales for gb1 vs rel-entropy estimate 
    for tau in imp_times_star:
        plt.plot([1,2],[tau,tau], 'k-')
    for tau in imp_times:
        plt.plot([3,4],[tau,tau], 'r-')
    plt.xlim(0,5)
    plt.xticks([1.5, 3.5], ['Fs-6383', identifier])
    plt.yscale('log')
    plt.ylabel('implied timescales (ns)')
    panel += 1


plt.tight_layout()
outfile = 'Fs-6383_vs_%s.pdf'%identifier
#print 'Writing', outfile, '(may take a minute)...'
#plt.show()
plt.savefig(outfile)
#print '...Done.'
