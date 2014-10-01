import numpy as np
import scipy
from scipy.optimize import minimize


def generate_W_sym(quantity,gamma):
    row=col=len(quantity)
    w = np.zeros((row,col))
    for i in range(row):
        for j in range(col):
            w[i,j] = np.abs(quantity[i] - quantity[j])
    w_sym = 0.5*(w+w.transpose())
    W_sym = np.exp(-gamma*w_sym)

    return W_sym,w

def solve_beta_lambda(W_sym,population,sigma=1e-5):
    beta= 1.0*np.ones((W_sym.shape[0],1))
    while True:
        lam = population.reshape(-1,1)/np.dot(W_sym,beta)
        beta_new = population.reshape(-1,1)/np.dot(W_sym,lam)
        lam_new = population.reshape(-1,1)/np.dot(W_sym,beta_new)
        if np.abs(beta_new-beta).max()<=sigma and np.abs(lam_new-lam).max()<=sigma:
            break
        else:
            beta = beta_new

    return beta_new,lam_new


def cal_aver_ensemble(populaiton,w,K):
    aver_w = 0
    for i in range(K.shape[0]):
        for j in range(K.shape[1]):
            aver_w = aver_w + populaiton[i]*K[i,j]*w[i,j]
    print aver_w

def estimate_K(W_sym,population,rho):
    K =np.zeros((W_sym.shape))
    for i in range(W_sym.shape[0]):
        for j in range(W_sym.shape[1]):
            K[i,j] = rho[i]*rho[j]*W_sym[i,j]/population[i]
    return K


p = np.loadtxt('Populations.dat')
d = np.loadtxt('mean_Nv_state.dat')
print p.shape

W_sym,w = generate_W_sym(d,1)
beta,lam = solve_beta_lambda(W_sym,p,sigma=1e-5)
rho = np.sqrt(beta*lam)
K = estimate_K(W_sym,p,rho)
aver_w = cal_aver_ensemble(W_sym,w,K)
print "K:",K
print "aver_w",aver_w




