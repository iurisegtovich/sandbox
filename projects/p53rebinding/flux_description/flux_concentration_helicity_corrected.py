from flux_utils import *
import numpy as np
import matplotlib.pyplot as plt

p = np.loadtxt("Populations.dat")

h_sim = p[2]/(p[2]+p[3])

print h_sim
print p[0]/p[0]+p[1]
h_exp = 0.64

c_p53 = 7.1*10**-3 #Mol
c_mdm2 = 7.1*10**-3

k_off_t = 1.8e+04
k_on_t = 8.7e+07
k_d_t = k_off_t/k_on_t

k_off_w = 3.2e+06
k_on_w = 2.1e+09
k_d_w = k_off_w/k_on_w

k_wt_l = 1.1e+06

k_wt = 1.8e+04
k_tw = 1.4e+06
k_eq_wt = k_wt/k_tw

k_off_overall = 8.6e+05 # 1/s
k_off_exp = 2.0

k_on_overall = 5.8e+07 #1/(M*s)
k_on_exp = 9.2e+06

factor_k_off =  k_off_exp/k_off_overall
factor_k_on = k_on_exp /k_on_overall

k_off_t = k_off_t*factor_k_off
k_on_t = k_on_t*factor_k_on
k_d_t = k_off_t/k_on_t
k_on_w = k_on_w*factor_k_on
k_off_w = k_off_w*factor_k_off
k_d_w = k_off_w/k_on_w



k_wt = k_wt*(h_exp*(1-h_sim)/(h_sim*(1-h_exp)))**0.5
k_tw = k_tw*(h_sim*(1-h_exp)/(h_exp*(1-h_sim)))**0.5

k_eq_wt = k_wt/k_tw

concentrations = 10**(np.arange(-7,0,0.01))
concentrations = np.insert(concentrations,0,7.1e-3)
concentrations = sorted(concentrations)

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
    plt.show()

fn1 = "flux_ratio_corrected_p53(mM)_vs_mdm2_helicity%d.txt"%(int(h_exp*100))
fn2 = "flux_ratio_corrected_p53(uM)_vs_mdm2_helicity%d.txt"%(int(h_exp*100))
fn3 = "binding_affinity_corrected_p53(mM)_vs_mdm2_helicity%d.txt"%(int(h_exp*100))
np.savetxt(fn1,flux_ratio_vs_mdm2)
np.savetxt(fn2,flux_ratio1_vs_mdm2)
np.savetxt(fn3,binding_affinity_vs_mdm2)
print "Wrote: %s,%s,%s"%(fn1,fn2,fn3)


