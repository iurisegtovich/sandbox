import os, sys, glob
import numpy as np

def estimate_perturbed_rates(T_star, pi_star, pi, tol=1e-12, maxsteps=50000, verbose=True):
    """Use the minimum relative-entropy estimator to estimate
    new rates upon perturbation of an MSM.
    
    NOTE:  So far this only works for *dense* matrices! --VAV
    
    INPUTS
    T_star    - a matrix of unperturbed transition probabilities $p_{ij}^*$
    pi_star   - a vector of unperturbed equilibirum populations $\pi_i^*$
    pi        - a vector of *perturbed* equilibirum populations $\pi_i$
    
    OPTIONS
    tol       - the desired tolerance to achieve numerical convergence (Default: 1e-12)
    maxsteps  - the maximum number of iterations to attempt
    
    RETURNS
    T         - a matrix of unperturbed transition probabilities p_ij
    """
    
    nstates = T_star.shape[0]

    print 'T_star', T_star
    
    # initial guess at w_i = exp(-v_i/pi_star)
    w = np.ones(nstates)

    T = np.zeros(T_star.shape)
    old_T = T.copy()

    for step in range(maxsteps):

        # estimate lagrange multipliers exp_negdelta_ij
        exp_negdelta_ij = np.zeros(T_star.shape)
        for i in range(nstates):
            for j in range(nstates):
                if T_star[i,j] != 0.:
                    exp_negdelta_ij[i,j] = np.sqrt( (pi[j]*T_star[j,i]*w[j]) / (pi[i]*T_star[i,j]*w[i])  )   

        # from the exp_negdelta_ij, estimate p_ij's
        for i in range(nstates):
            for j in range(nstates):
                T[i,j] = T_star[i,j]*exp_negdelta_ij[i,j]*w[i]

        # print 'T' , T
        # print 'old_T', old_T
        
        # calculate the convergence error
        err = np.max( np.abs(T - old_T))
        # print 'err', err 
        if err < tol:
            print 'Reached convergence tolerance |T - oldT| < %e'%tol
            break

        if verbose:
            print 'step', step, 'err', err


        # use these value of exp_negdelta_ij to update weights w[i] 
        for i in range(nstates):
            w[i] = w[i]/(T[i,:]).sum() 

        old_T = T.copy()  # keep track of last iteration so we can compare

    return T
    


