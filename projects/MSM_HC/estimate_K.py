import numpy as np


def generate_W_sym(quantity,gamma):
    row=col=len(quantity)
    W = np.zeros((row,col))
    for i in range(row):
        for j in range(col):
            wij = np.abs(quantity[i] - quantity[j])
            W[i,j] = np.exp(-gamma*wij)
    return np.matrix(W)

def solve_rho(W_sym,population,sigma=0.001):

    rho=np.ones((W_sym.shape[0],1))
    while True:
        rho_new = W_sym.I*(population.reshape(-1,1)/rho)
        print rho_new[0]
        if np.abs(rho_new-rho).max()<=sigma:
            break
        else:
            rho = rho_new
    return rho_new






p = np.loadtxt('Populations.dat')
d = np.loadtxt('mean_Nv_state.dat')


W_sym = generate_W_sym(d,1)
rho = solve_rho(W_sym,p,sigma=0.001)
print rho





