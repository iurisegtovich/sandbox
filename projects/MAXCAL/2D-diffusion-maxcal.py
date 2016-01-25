import numpy as np
import matplotlib.pyplot as plt

pid = 6383

T = np.load("./data/tProb-%d.npy"%pid)
C_sym = np.load("./data/tCounts-%d.npy"%pid)
m,n = T.shape
obs = np.loadtxt("./data/Nhs_state-%d.dat"%pid)

# Equilibrium populations should be uniform
if (0):
    pi_a = np.ones(T.shape[0])/float(nstates)
else:
    ### The stationary eigenvector of np.transpose(T) seems to work a bit better
    evals, evecs = np.linalg.eig(np.transpose(T))
    pi_a = evecs[:,np.argmax(evals)]
    pi_a = pi_a/pi_a.sum()
    print 'pi_a.sum()', pi_a.sum(), 'pi_a', pi_a


# Calculate matrix of jump indicators, N_ab
N_ab = np.ones(T.shape)
for i in range(T.shape[0]):
    N_ab[i,i] = 0.

# ... and the mean jump rate from the observed counts
N_mean = (N_ab*C_sym).sum()/C_sym.sum()
print 'The mean jump rate N_mean =', N_mean


# Calculate the matrix of dynamical observables r_ab
### in our case: the distance between state centers
r_ab = np.zeros(T.shape)
nstates = len(obs)
for a in range(nstates):
    for b in range(nstates):
        r_ab[a,b] = np.abs(obs[a]-obs[b])

# ... and the mean r_ab from the observed counts
r_mean = (r_ab*C_sym).sum()/C_sym.sum()
print 'The mean distance traveled per unit time, r_mean =', r_mean


plt.figure(figsize=(10,4))
plt.subplot(1,2,1)
plt.pcolor(N_ab)

plt.subplot(1,2,2)
plt.pcolor(r_ab)
plt.show()

def D(x, pi):
    """$\mathcal{D}$ is the non-linear operator $\mathcal{D}[\vec{x}]_a = \pi_a/x_a$"""
    return pi * (1./x)

##############
### Let's do a scan for the best alpha and gamma

alpha_values = np.arange(-2.5, -0.1,0.01)
gamma_values = np.arange(0.4, 2.8,0.01)

# store the squared error (N_maxcal - N_mean)**2 + (r_maxcal - r_mean)**2 for each parameter set
chi2_results = np.zeros( (alpha_values.shape[0],gamma_values.shape[0]))

besti, bestj, best_chi2 = 0,0, 1.0e99

for i in range(len(alpha_values)):
    for j in range(len(gamma_values)):


        alpha = alpha_values[i]
        gamma = gamma_values[j]


        # calculate W_ab
        W = np.exp(-alpha*N_ab) * np.exp(-gamma*r_ab)

        tol = 1.0e-12
        max_delta = 1.0
        trial = 0
        beta_a = np.ones(nstates) # initial guess
        lambda_b = np.ones(nstates) # initial guess
        while max_delta > tol:

            new_beta_a = D( np.dot(W,lambda_b), pi_a)
            new_lambda_b = D( np.dot(np.transpose(W), new_beta_a), pi_a)

            # measure the difference between old and new estimates
            diff_beta = new_beta_a - beta_a
            delta_beta = np.sqrt(np.dot(diff_beta,diff_beta))

            diff_lambda = new_lambda_b - lambda_b
            delta_lambda = np.sqrt(np.dot(diff_lambda,diff_lambda))

            max_delta = max(delta_beta,delta_lambda)
            #print 'trial', trial, 'beta_a =', beta_a[0:3], '... lambda_b =', lambda_b[0:3],'...', 'delta_beta', delta_beta, 'delta_lambda', delta_lambda

            beta_a = new_beta_a
            lambda_b = new_lambda_b
            trial += 1

        # print 'trial', trial, 'beta_a =', beta_a[0:3], '... lambda_b =', lambda_b[0:3],'...', 'delta_beta', delta_beta, 'delta_lambda', delta_lambda

        # Use the converged values of beta_a and lambda_b to get T^{maxcal}_ab
        T_maxcal = np.zeros(T.shape)
        for a in range(nstates):
            for b in range(nstates):
                T_maxcal[a,b] = beta_a[a]/pi_a[a]*lambda_b[b]*W[a,b]
        # print 'T_maxcal', T_maxcal

        # print '#### For alpha =', alpha, 'gamma =', gamma, '####'
        N_maxcal = np.dot(pi_a, (T_maxcal*N_ab).sum(axis=1))
        # print 'N_maxcal =',N_maxcal, 'N_mean = ', N_mean
        r_maxcal = np.dot(pi_a, (T_maxcal*r_ab).sum(axis=1))
        #print 'r_maxcal =',r_maxcal, 'r_mean = ', r_mean

        chi2 = (N_maxcal - N_mean)**2 + (r_maxcal - r_mean)**2
        #print 'chi2 =', chi2

        if chi2 < best_chi2:
            best_chi2 = chi2
            besti, bestj = i,j
            best_alpha = alpha_values[besti]
            best_gamma = gamma_values[bestj]
            best_T_maxcal = T_maxcal
            print '#### BEST so far: alpha =', alpha, 'gamma =', gamma, '####'
            print 'N_maxcal =',N_maxcal, 'N_mean = ', N_mean
            print 'r_maxcal =',r_maxcal, 'r_mean = ', r_mean
            print 'chi2 =', chi2


        chi2_results[i,j] = chi2

plt.figure( figsize=(6,6))
X,Y = np.meshgrid(alpha_values, gamma_values)
plt.contour(X,Y, np.transpose(np.log(chi2_results)))
plt.xlabel('$\\alpha$')
plt.ylabel('$\\gamma$')
plt.show()
