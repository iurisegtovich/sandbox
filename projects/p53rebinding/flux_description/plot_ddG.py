import numpy as np
import matplotlib.pyplot as plt

concentrations = 10**(np.arange(-7,0,0.01))
concentrations = np.insert(concentrations,0,7.1e-3)
concentrations = sorted(concentrations)

ddG_28 = np.loadtxt("ddG_28_h_corrected_rates.txt")
ddG_64 = np.loadtxt("ddG_64_h_corrected_rates.txt")

for i,c in enumerate(concentrations):
    if abs(c-7.1e-3)<1e-6:
        print ddG_28[i]
	print ddG_64[i]
	print "difference:",ddG_64[i]-ddG_28[i]

plt.figure(figsize=(8,6))
plt.plot(concentrations,ddG_28)
plt.plot(concentrations,ddG_64)
plt.legend(["p53=7.1 mM h%=28%","p53=7.1 mM h%=64%"],loc='best')
params = {'legend.fontsize': 15,
          'legend.handlelength': 2}
plt.rcParams.update(params)
plt.tick_params(axis='both', which='major', labelsize=15)
plt.xlabel("MDM2 concentration(M)",fontsize=20)
plt.ylabel(r"$\Delta \Delta$ G (kcal/mol)",fontsize=20)
plt.xscale('log')
plt.plot([0.0071, 0.0071], [-1.0,-4.5], 'k--')
plt.ylim([-4.5,-1.0])
plt.margins(0.10)
fig_fn = "ddG_corrected_rates.pdf"
plt.savefig(fig_fn)
print "Saved: %s"%fig_fn
