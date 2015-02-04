import numpy as np
import scipy
from scipy.optimize import minimize
from scipy.io import mmwrite,mmread


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


def cal_aver_ensemble(population,w,K):
    aver_w = 0
    for i in range(K.shape[0]):
        for j in range(K.shape[1]):
            aver_w = aver_w + population[i]*K[i,j]*w[i,j]
    return aver_w

def estimate_K(W_sym,population,rho):
    K =np.zeros((W_sym.shape))
    for i in range(W_sym.shape[0]):
        for j in range(W_sym.shape[1]):
            K[i,j] = rho[i]*rho[j]*W_sym[i,j]/population[i]
    return K


p = np.loadtxt('Populations.LifsonRoig.dat')
d = np.loadtxt('mean_Nv_state.dat')
"""
aver_w_all = []
gamma = 0.05
while gamma <=0.2:
    W_sym,w = generate_W_sym(d,gamma)
    print gamma
    #print "W_sym",W_sym
    #print "w",w
    beta,lam = solve_beta_lambda(W_sym,p,sigma=1e-5)
    rho = np.sqrt(beta*lam)
    K = estimate_K(W_sym,p,rho)
    aver_w = cal_aver_ensemble(p,w,K)
    aver_w_all.append(aver_w)
    print aver_w
    #print "K:",K
    #print "aver_w",aver_w
    #mmwrite('K_lifsonroig.gamma%d.mtx'%gamma,K)
    gamma = gamma + 0.01
np.savetxt('aver_w.smallgamma.dat',aver_w_all)
"""
#aver_w from MD
#aver_w 2.23330143749
#std 1.79336175161
# aver_w from MaxCal
#gamma = 0.12
#aver_w = 2.23220566646
W_sym,w = generate_W_sym(d,0.12)
#print "W_sym",W_sym
#print "w",w
beta,lam = solve_beta_lambda(W_sym,p,sigma=1e-5)
rho = np.sqrt(beta*lam)
K = estimate_K(W_sym,p,rho)
aver_w = cal_aver_ensemble(p,w,K)
print aver_w
mmwrite('K_lifsonroig.gamma_0.12.mtx',K)




