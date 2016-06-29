import numpy as np
import matplotlib.pyplot as plt

flux_raw_data = np.loadtxt("flux_ratio_corrected_p53(mM)_vs_mdm2.txt")
#flux_raw_data_p53um = np.loadtxt("flux_ratio_corrected_p53(uM)_vs_mdm2.txt")
flux_helicity28 = np.loadtxt("flux_ratio_corrected_p53(mM)_vs_mdm2_helicity28.txt")
#flux_helicity28_p53um = np.loadtxt("flux_ratio_corrected_p53(uM)_vs_mdm2_helicity28.txt")
flux_helicity64 = np.loadtxt("flux_ratio_corrected_p53(mM)_vs_mdm2_helicity64.txt")
#flux_helicity64_p53um = np.loadtxt("flux_ratio_corrected_p53(uM)_vs_mdm2_helicity64.txt")

concentrations = 10**(np.arange(-7,0,0.01))
concentrations = np.insert(concentrations,0,7.1e-3)
concentrations = sorted(concentrations)

if 1:
    plt.figure(figsize=(8, 6))
#    plt.plot(concentrations,flux_raw_data_p53um)
    plt.plot(concentrations,flux_raw_data)
#    plt.plot(concentrations,flux_helicity28_p53um)
    plt.plot(concentrations,flux_helicity28)
#    plt.plot(concentrations,flux_helicity64_p53um)
    plt.plot(concentrations,flux_helicity64)
    plt.ylabel(r"$\frac{F_{cs}}{F_{cs}+F_{if}}$",fontsize=25)
    plt.xlabel(r"MDM2 concentration (M)",fontsize=20)
    params = {'legend.fontsize': 12,
          'legend.handlelength': 2}
    plt.rcParams.update(params)
    plt.tick_params(axis='both', which='major', labelsize=15)
#    plt.legend(["p53=7.1uM h%=0.11%","p53=7.1mM h%=0.11%","p53=7.1uM h%=28%","p53=7.1mM h%=28%",
#                "p53=7.1uM h%=64%","p53=7.1mM h%=64%"], loc='best')
    plt.legend(["p53=7.1mM h%=0.11%","p53=7.1mM h%=28%",
                "p53=7.1mM h%=64%"], loc='best')
    plt.plot([0.0071, 0.0071], [0, 0.5], 'k--')
    plt.xscale('log')
    plt.ylim([0,0.5])
    plt.margins(0.10)
    plt.subplots_adjust(left=0.20)
    figfn = "figure_flux_helicity_corrected.pdf"
    plt.savefig(figfn)
    print "Saved: %s"%figfn
