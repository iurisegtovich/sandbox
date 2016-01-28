import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import distance
import scipy.linalg
from msmbuilder.msm.core import _solve_msm_eigensystem

pid = 6383

T = np.load("./data/tProb-%d.npy"%pid)
C_sym = np.load("./data/tCounts-%d.npy"%pid)
m,n = T.shape
obs = np.loadtxt("./data/Nhs_state-%d.dat"%pid)
rmsd = np.load("./data/RMSDs-states-macro40-p%d.npy"%pid)
tica = np.load("./data/tica_8-states-macro40-p%d.npy"%pid)

# Equilibrium populations should be uniform
if (0):
    pi_a = np.ones(T.shape[0])/float(nstates)
else:
    ### The stationary eigenvector of np.transpose(T) seems to work a bit better
    evals, evecs = np.linalg.eig(np.transpose(T))
    pi_a = evecs[:,np.argmax(evals)]
    pi_a = pi_a/pi_a.sum()
    print 'pi_a.sum()', pi_a.sum(), 'pi_a', pi_a
    lagtime=50
    impliedts = []
    top_ten_evals = evals[np.argsort(-evals)[1:11]]
    for i in top_ten_evals:
        impliedts.append(-lagtime/np.log(i))
    print "Implied Timescales(MSM):",impliedts


# Calculate matrix of jump indicators, N_ab
N_ab = np.ones(T.shape)
for i in range(T.shape[0]):
    for j in range(T.shape[1]):
        if T[i,j] == 0:
            N_ab[i,j] = 0
for i in range(T.shape[0]):
    N_ab[i,i] = 0.

# ... and the mean jump rate from the observed counts
N_mean = (N_ab*C_sym).sum()/C_sym.sum()
print 'The mean jump rate N_mean =', N_mean


# Calculate the matrix of dynamical observables r_ab
### in our case: the distance between state centers
r_ab = np.zeros(T.shape)
nstates = len(tica)
for a in range(nstates):
    for b in range(nstates):
        #r_ab[a,b] = np.abs(rmsd[a]-rmsd[b])
        r_ab[a,b] = distance.euclidean(tica[a][:],tica[b][:])

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
#for rmsd
#alpha_values = np.arange(0.81, 0.825,0.001)
#gamma_values = np.arange(9.75, 9.8,0.001)

#for 8-tics
alpha_values = np.arange(0.6, 0.7,0.005)
gamma_values = np.arange(0.88, 0.90,0.005)

#for 2-tics
#alpha_values = np.arange(1.0, 1.6,0.1)
#gamma_values = np.arange(1.5, 3.0,0.1)

# store the squared error (N_maxcal - N_mean)**2 + (r_maxcal - r_mean)**2 for each parameter set
chi2_results = np.zeros( (alpha_values.shape[0],gamma_values.shape[0]))

besti, bestj, best_chi2 = 0,0, 1.0e99

for i in range(len(alpha_values)):
    for j in range(len(gamma_values)):


        alpha = alpha_values[i]
        gamma = gamma_values[j]


        # calculate W_ab
        W = np.exp(-alpha*N_ab) * np.exp(-gamma*r_ab)

        tol = 1.0e-8
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
            #print "delta_beta",delta_beta,"delta_lambda",delta_lambda
            #print "max_delta",max_delta
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

plt.figure(figsize=(6,6))
for a in range(T.shape[0]):
    for b in range(T.shape[1]):
        plt.plot(T[a,b], best_T_maxcal[a,b], '.' )
plt.xlabel('$T_{ab}$ from MSM', fontsize=12)
plt.ylabel('maximum caliber $T_{ab}$', fontsize=12)
plt.xscale('log')
plt.yscale('log')
plt.ylim([1e-8,1])
plt.xlim([1e-8,1])
#plt.plot([1,1e-3],[1,1e-3],'k-')
plt.plot([1,1e-8],[1,1e-8],'k-')
#figfn="rmsd-maxcal.pdf"
figfn = "tICs-8-maxcal.pdf"
plt.title(figfn)
plt.savefig(figfn)
plt.show()

#plot colormap of T and T_maxcal
plt.figure(figsize=(12,5))
plt.subplot(1,2,1)
plt.pcolor(T)
plt.title("Color map of T_MSM")

plt.subplot(1,2,2)
plt.pcolor(best_T_maxcal)
plt.title("Color map of T_MaxCal")

plt.savefig("pcolor_T_8_tICs.png")
plt.show()

#calculate implied timescales associated with T_maxcal
u_maxcal, lv_maxcal, rv_maxcal = scipy.linalg.eig(np.transpose(best_T_maxcal),left=True,right=True)
print lv_maxcal[:,np.argmax(u_maxcal)]
print rv_maxcal[:,np.argmax(u_maxcal)]
print u_maxcal[np.argmax(u_maxcal)]
ind_maxcal = np.argsort(-u_maxcal)
top_ten_evals = u_maxcal[ind_maxcal[1:11]]

lagtime = 50
impliedts = []
for i in top_ten_evals:
    impliedts.append(-lagtime/np.log(i))
print "Implied Timescales(MaxCal):",impliedts

u_msm, lv_msm, rv_msm = scipy.linalg.eig(np.transpose(T),left=True,right=True)
print lv_msm[:,np.argmax(u_msm)]
print rv_msm[:,np.argmax(u_msm)]
impliedts = []
ind_msm = np.argsort(-u_msm)
top_ten_evals = u_msm[ind_msm[1:11]]
for i in top_ten_evals:
    impliedts.append(-lagtime/np.log(i))
print "Implied Timescales(MSM):",impliedts

p = 1./40.*np.ones(40)
"""
plt.figure(figsize=(12,12))

for j in range(5):
    ax = plt.subplot(5,2,j*2+1)
    #plt.plot(range(len(msm.left_eigenvectors_[:,j])),msm.left_eigenvectors_[:,j],"ro--")
    #plt.plot(range(len(msm.left_eigenvectors_[:,j])),[0]*len(msm.left_eigenvectors_[:,j]),"k--")
    plt.hold(True)

    if np.dot(p,lv_msm[:,j]) > 0 :
        for stateid in range(40):
            plt.vlines(stateid,min(0,rv_msm[stateid,j]),max(rv_msm[stateid,j],0),linestyles='solid',linewidth="2")
    else:
        for stateid in range(40):
            plt.vlines(stateid,min(0,-rv_msm[stateid,j]),max(-rv_msm[stateid,j],0),linestyles='solid',linewidth="2")

    if j == 0:
        plt.title("Fs-EEE-msm")

for j in range(5):
    ax = plt.subplot(5,2,j*2+2)
    #plt.plot(range(len(msm.left_eigenvectors_[:,j])),msm.left_eigenvectors_[:,j],"ro--")
    #plt.plot(range(len(msm.left_eigenvectors_[:,j])),[0]*len(msm.left_eigenvectors_[:,j]),"k--")
    plt.hold(True)
    if np.dot(p,lv_maxcal[:,j]) > 0 :
        for stateid in range(40):
            plt.vlines(stateid,min(0,rv_maxcal[stateid,j]),max(rv_maxcal[stateid,j],0),linestyles='solid',linewidth="2")
    else:
        for stateid in range(40):
            plt.vlines(stateid,min(0,-rv_maxcal[stateid,j]),max(-rv_maxcal[stateid,j],0),linestyles='solid',linewidth="2")
    if j == 0:
        plt.title("Fs-EEE-maxcal")
plt.savefig("eigenvectors-msm-maxcal.pdf")
plt.show()
"""

k = 11

u, lv, rv = _solve_msm_eigensystem(T, k)
V = rv
S = np.diag(lv[:,0])
C = S.dot(T)
try:
    trace = np.trace(V.T.dot(C.dot(V)).dot(np.linalg.inv(V.T.dot(S.dot(V)))))
except np.linalg.LinAlgError:
    trace = np.nan
print "T_MSM score",trace

u, lv, rv = _solve_msm_eigensystem(best_T_maxcal, k)
V = rv
S = np.diag(lv[:,0])
C = S.dot(T)
try:
    trace = np.trace(V.T.dot(C.dot(V)).dot(np.linalg.inv(V.T.dot(S.dot(V)))))
except np.linalg.LinAlgError:
    trace = np.nan
print "T_Maxcal score",trace

