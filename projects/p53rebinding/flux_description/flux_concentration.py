from flux_utils import *
import numpy as np
import matplotlib.pyplot as plt


c_p53 = 7.1*10**-3 #Mol
c_mdm2 = 7.1*10**-3

k_off_t = 1.8e+04
k_on_t = 8.7e+07
#k_on_t = 6.0e+11
k_d_t = k_off_t/k_on_t

k_off_w = 3.2e+06
k_on_w = 2.1e+09
#k_on_w = 1.2e+10
k_d_w = k_off_w/k_on_w

k_wt_l = 1.1e+06

k_wt = 1.8e+04
k_tw = 1.4e+06
k_eq_wt = k_wt/k_tw

concentrations = 10**(np.arange(-7,0,0.1))

flux_ratio_vs_mdm2 = []
binding_affinity_vs_mdm2 = []

for c_mdm2 in concentrations:
    mdm2_free = solve_free_protein(mdm2_tot=c_mdm2, p53_tot=7.1e-3,
                         k_d_w=k_d_w, k_d_t=k_d_t,
                         k_eq_wt=k_eq_wt)
    flux_cs,flux_if = estimate_flux(p53_tot=7.1e-3,k_wt=k_wt
                    ,k_d_w=k_d_w,k_d_t=k_d_t
                    ,k_on_t=k_on_t,k_on_w=k_on_w
                    ,k_wt_l=k_wt_l,k_eq_wt=k_eq_wt,mdm2_free=mdm2_free)
    affinity = estimate_binding_affinity(p53_tot=7.1e-3,k_wt=k_wt,
                                         k_d_w=k_d_w,k_d_t=k_d_t,
                                         k_on_t=k_on_t,k_on_w=k_on_w,
                                         k_wt_l=k_wt_l,mdm2_free=mdm2_free,
                                         k_eq_wt=k_eq_wt)
    #print "F_cs:",flux_cs,"F_if:",flux_if,"F_cs/(F_cs+F_if):",flux_cs/(flux_cs+flux_if)
    flux_ratio_vs_mdm2.append(flux_cs/(flux_cs+flux_if))
    binding_affinity_vs_mdm2.append(affinity)

flux_ratio1_vs_mdm2 = []

for c_mdm2 in concentrations:
    mdm2_free = solve_free_protein(mdm2_tot=c_mdm2, p53_tot=7.1e-6,
                         k_d_w=k_d_w, k_d_t=k_d_t,
                         k_eq_wt=k_eq_wt)
    flux_cs,flux_if = estimate_flux(p53_tot=7.1e-6,k_wt=k_wt
                    ,k_d_w=k_d_w,k_d_t=k_d_t
                    ,k_on_t=k_on_t,k_on_w=k_on_w
                    ,k_wt_l=k_wt_l,k_eq_wt=k_eq_wt,mdm2_free=mdm2_free)
    #print "F_cs:",flux_cs,"F_if:",flux_if,"F_cs/(F_cs+F_if):",flux_cs/(flux_cs+flux_if)
    flux_ratio1_vs_mdm2.append(flux_cs/(flux_cs+flux_if))

if 0:
    plt.plot(concentrations,flux_ratio_vs_mdm2)
    plt.plot(concentrations,flux_ratio1_vs_mdm2)
    plt.ylabel(r"$\frac{F_{cs}}{F_{cs}+F_{if}}$",fontsize=25)
    plt.xlabel(r"MDM2 concentration (M)",fontsize=15)
    params = {'legend.fontsize': 12,
          'legend.handlelength': 2}
    plt.rcParams.update(params)
    plt.legend(["p53=7.1 mM","p53=7.1 uM"],loc='best')
    plt.xscale('log')
    #plt.savefig("flux_rmsd_2.0_population_weighted.pdf")
    plt.savefig("flux_rmsd_2.0.pdf")

