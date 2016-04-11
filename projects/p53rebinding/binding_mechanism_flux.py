import numpy as np
from scipy.optimize import broyden1,broyden2,newton
import matplotlib.pyplot as plt



def solve_free_protein(mdm2_tot,p53_tot,k_d_w,k_d_t,k_eq_wt,alpha=1.):

    # k_off_w/k_on_w = k_unfold_bound_to_unfold_unbound/k_unfold_unbound_to_unfold_bound
    #                = 7.120e+06/1.537e+07
    k_d_w *= alpha
    # k_off_t/k_on_t = k_fold_bound_to_fold_unbound/k_fold_unbound_to_fold_bound
    #                = 1.704e+06/6.565e+06
    k_d_t *= alpha
    # k_w_to_t/k_t_to_w = k_unfold_unbound_to_fold_unbound/k_fold_unbound_to_unfold_unbound
    #                = 6.145e+05/2.162e+06
    k_eq_wt *= alpha

    def F(x):
        return mdm2_tot*(1+x/k_d_w+k_eq_wt*(1+x/k_d_t))-p53_tot*(1/k_d_w+k_eq_wt/k_d_t)-x*(1+x/k_d_w+k_eq_wt*(1+x/k_d_t))

    #mdm2_free = broyden1(F,0.5*p53_tot)
    a = 1/k_d_w+k_eq_wt/k_d_t
    b = 1+k_eq_wt-mdm2_tot/k_d_w-k_eq_wt/k_d_t*mdm2_tot
    c = p53_tot/k_d_w+p53_tot*k_eq_wt/k_d_t-(1+k_eq_wt)*mdm2_tot
    print a,b,c
    mdm2_free = -b+(b*b-4*a*c)**0.5/(2*a)
    return mdm2_free

def estimate_flux(p53_w_unbound,k_eq_wt,k_d_w,k_d_t,mdm2_free):

    p53_t_unbound = k_eq_wt*p53_w_unbound
    p53_w_bound = p53_w_unbound/k_d_w*mdm2_free
    #p53_t_bound = k_eq_wt*p53_w_unbound/k_d_t*mdm2_free

    flux_cs = (1/(k_wt*p53_w_unbound)+1/(k_on_t*p53_t_unbound*mdm2_free))**-1
    flux_if = (1/(k_on_w*p53_w_unbound*mdm2_free+1/k_wt_bound*p53_w_bound))**-1

    return flux_cs,flux_if




c_p53 = 7.1*10**-3 #Mol

# k_off_w/k_on_w = k_unfold_bound_to_unfold_unbound/k_unfold_unbound_to_unfold_bound
#                = 7.120e+06/1.537e+07
k_d_w = 7.120e+06/8.980e+09
# k_off_t/k_on_t = k_fold_bound_to_fold_unbound/k_fold_unbound_to_fold_bound
#                = 1.704e+06/6.565e+06
k_d_t = 1.704e+06/4.709e+10
# k_w_to_t/k_t_to_w = k_unfold_unbound_to_fold_unbound/k_fold_unbound_to_unfold_unbound
#                = 6.145e+05/2.162e+06
k_eq_wt = 6.145e+05/2.162e+06

print solve_free_protein(mdm2_tot=c_p53, p53_tot=c_p53,
                         k_d_w=k_d_w, k_d_t=k_d_t,
                         k_eq_wt=k_eq_wt, alpha=1)


