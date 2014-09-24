import numpy as np


def generate_W_sym(quantity,gamma):
    row=col=len(quantity)
    W = np.zeros((row,col))
    for i in range(row):
        for j in range(col):
            wij = np.abs(quantity[i] - quantity[j])
            W[i,j] = np.exp(-gamma*wij)
    return np.matrix(W)

def solve_rho(W_sym,sigma=0.001):
    rho=np.ones(W_sym.shape[0])



p = np.loadtxt('Populations.dat')
d = np.loadtxt('mean_Nv_state.dat')


W = generate_W_sym(d,1)
print W
