import numpy as np
from scipy.optimize import broyden1,broyden2,newton
import matplotlib.pyplot as plt



def solve_free_protein(mdm2_tot,p53_tot,alpha=1.):

    # k_off_w/k_on_w = k_unfold_bound_to_unfold_unbound/k_unfold_unbound_to_unfold_bound
    #                = 7.120e+06/1.537e+07
    k_d_w = 7.120e+06/1.537e+07*alpha
    # k_off_t/k_on_t = k_fold_bound_to_fold_unbound/k_fold_unbound_to_fold_bound
    #                = 1.704e+06/6.565e+06
    k_d_t = 1.704e+06/6.565e+06*alpha
    # k_w_to_t/k_t_to_w = k_unfold_unbound_to_fold_unbound/k_fold_unbound_to_unfold_unbound
    #                = 6.145e+05/2.162e+06
    k_eq_wt = 6.145e+05/2.162e+06*alpha

    def F(x):
        return mdm2_tot*(1+x/k_d_w+k_eq_wt*(1+x/k_d_t))-(p53_tot)*(1/k_d_w+k_eq_wt/k_d_t)-x*(1+x/k_d_w+k_eq_wt*(1+x/k_d_t))

    mdm2_free = newton(F,0.1*p53_tot)

    return mdm2_free

def estimate_flux(c_free_protein):
    


c_p53 = 7.1*10**-3 #Mol

print solve_free_protein(mdm2_tot=c_p53,p53_tot=c_p53,alpha=0.01)

