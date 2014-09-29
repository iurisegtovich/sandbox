import numpy as np
import scipy
from scipy.optimize import minimize


def generate_W_sym(quantity,gamma):
    row=col=len(quantity)
    W = np.zeros((row,col))
    for i in range(row):
        for j in range(col):
            wij = np.abs(quantity[i] - quantity[j])
            W[i,j] = np.exp(-gamma*wij)

    return W

def solve_rho(W_sym,population,sigma=0.001):
    # self-consistently slove the equation. Doesn't work yet.
    rho=5.0*np.ones(W_sym.shape[0])
    i=0
    while i<100:
        i = i+1
        v = np.dot(rho,W_sym)
        rho_new = population/v
        print "step %d"%i,"rho:",rho
        if np.abs(rho_new-rho).max()<=sigma:
            break
        else:
            rho = rho_new
    return rho_new

def solve_rho_min(W_sym,population):
    
    fun = lambda rho: np.sum((rho*np.dot(W_sym,rho) - population.reshape(-1,1))**2)
    #rho= 1.0*np.ones((W_sym.shape[0],1))
    rho = np.sqrt(population.reshape(-1,1))
    cons = ({'type':'ineq','fun':lambda rho: rho})
    bnds = tuple(bnds)
    res = minimize(fun,rho,method='SLSQP',bounds=bnds,options={'ftol': 1e-6, 'eps': 1e-6})
    print res
    return res['x']

def cal_aver_ensemble(W_sym,quantity,rho):
    aver_w = 0
    for i in range(W_sym.shape[0]):
        for j in range(W_sym.shape[1]):
            aver_w = aver_w + rho[i]*rho[j]*W_sym[i,j]*np.abs(quantity[i]-quantity[j])
    print aver_w

def F(rho):
    return rho*np.dot(W_sym,rho) - p.reshape(-1,1)



global p
p = np.loadtxt('recover/populations.dat')
#d = np.loadtxt('mean_Nv_state.dat')
d = np.arange(0,11)
print p.shape

global W_sym
W_sym = generate_W_sym(d,1)
#rho = solve_rho(W_sym,p)

rho_init= 1.0*np.ones((W_sym.shape[0],1))
print rho_init*np.dot(W_sym,rho_init) - p.reshape(-1,1)
rho = scipy.optimize.broyden1(F,rho_init,f_tol=1e-14)
print 'fun:',rho*np.dot(W_sym,rho) - p.reshape(-1,1)
print rho
#cal_aver_ensemble(W_sym,d,rho)



