import matplotlib
matplotlib.use("Agg")
from pylab import *
from scipy import stats

# Plot parameters
params = {
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "font.size": 9,
    "legend.fontsize": 9,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "text.usetex": True,
    "figure.figsize": (3.15, 3.15),
    "figure.subplot.left": 0.15,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.13,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0.15,
    "figure.subplot.hspace": 0.12,
    "lines.markersize": 6,
    "lines.linewidth": 2.0,
    "text.latex.unicode": True,
}
rcParams.update(params)
rc("font", **{"family": "sans-serif", "sans-serif": ["Times"]})

# Equations of state
eos_SF_rho = np.logspace(-10, 5, 1000)
eos_SF_T = (eos_SF_rho / 10**(-1))**(1./3.) * 8000.

# Plot the phase space diagram
figure()
subplot(111, xscale="log", yscale="log")
plot(eos_SF_rho, eos_SF_T, 'k-', lw=1.)

plot([1e-10, 1e-1], [8000, 8000], 'k:', lw=0.6)
plot([1e-1, 1e-1], [20, 8000], 'k:', lw=0.6)
plot([1e-1, 1e1], [20000., 20000.*10.**(2./3.)], 'k--', lw=0.6)
text(1e-1, 200000, "$n_{\\rm H}~\\widehat{}~{\\tt EOS\\_gamma\\_effective}$", va="top", rotation=43, fontsize=7)
text(0.95e-1, 25, "${\\tt EOS\\_density\\_norm\\_H\\_p\\_cm3}$", rotation=90, va="bottom", ha="right", fontsize=7)
text(5e-8, 8400, "${\\tt EOS\\_temperature\\_norm\\_K}$", va="bottom", fontsize=7)

scatter([1e-1], [8000], s=4, color='k')

xlabel("${\\rm Hydrogen~number~density}~n_{\\rm H}~[{\\rm cm^{-3}}]$", labelpad=0)
ylabel("${\\rm Temperature}~T~[{\\rm K}]$", labelpad=2)
xlim(3e-8, 3e3)
ylim(20., 2e5)
savefig("EAGLE_SF_EOS.png", dpi=200)
