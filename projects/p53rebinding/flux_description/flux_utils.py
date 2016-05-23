

def solve_free_protein(mdm2_tot,p53_tot,k_d_w,k_d_t,k_eq_wt):

    a = 1./k_d_w+k_eq_wt/k_d_t
    b = 1.+k_eq_wt-mdm2_tot/k_d_w-k_eq_wt/k_d_t*mdm2_tot+p53_tot/k_d_w+p53_tot*k_eq_wt/k_d_t
    c = -(1.+k_eq_wt)*mdm2_tot
    mdm2_free1 = (-b+(b*b-4.*a*c)**0.5)/(2.*a)
    mdm2_free2 = (-b-(b*b-4.*a*c)**0.5)/(2.*a)

    return max(mdm2_free1,mdm2_free2)

def estimate_flux(p53_tot,k_wt,k_d_w,k_d_t,k_on_t,k_on_w,k_wt_l,k_eq_wt,mdm2_free):

    p53_w_unbound = p53_tot/(1+mdm2_free/k_d_w+k_eq_wt+k_eq_wt/k_d_t*mdm2_free)
    p53_t_unbound = k_eq_wt*p53_w_unbound
    p53_w_bound = p53_w_unbound/k_d_w*mdm2_free
    #p53_t_bound = k_eq_wt*p53_w_unbound/k_d_t*mdm2_free

    flux_cs = (1./(k_wt*p53_w_unbound)+1./(k_on_t*p53_t_unbound*mdm2_free))**-1
    flux_if = (1./(k_on_w*p53_w_unbound*mdm2_free)+1./(k_wt_l*p53_w_bound))**-1

    return flux_cs,flux_if

def estimate_binding_affinity(p53_tot,k_wt,k_d_w,k_d_t,k_on_t,k_on_w,k_wt_l,mdm2_free,k_eq_wt):

    p53_w_unbound = p53_tot/(1+mdm2_free/k_d_w+k_eq_wt+k_eq_wt/k_d_t*mdm2_free)
    p53_t_unbound = k_eq_wt*p53_w_unbound
    p53_w_bound = p53_w_unbound/k_d_w*mdm2_free
    p53_t_bound = k_eq_wt*p53_w_unbound/k_d_t*mdm2_free
    
    return p53_t_bound/(p53_w_unbound+p53_t_unbound+p53_w_bound)



