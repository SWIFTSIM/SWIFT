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
eos_cool_rho = np.logspace(-5, 5, 1000)
eos_cool_T = eos_cool_rho**0. * 8000.
eos_Jeans_rho = np.logspace(-1, 5, 1000)
eos_Jeans_T = (eos_Jeans_rho/ 10**(-1))**(1./3.) * 4000.

# Plot the phase space diagram
figure()
subplot(111, xscale="log", yscale="log")
plot(eos_cool_rho, eos_cool_T, 'k-', lw=1.)
plot(eos_Jeans_rho, eos_Jeans_T, 'k-', lw=1.)
plot([1e-10, 1e-5], [8000, 8000], 'k:', lw=0.6)
plot([1e-10, 1e-1], [4000, 4000], 'k:', lw=0.6)
plot([1e-5, 1e-5], [20, 8000], 'k:', lw=0.6)
plot([1e-1, 1e-1], [20, 4000], 'k:', lw=0.6)
plot([3e-6, 3e-4], [28000, 28000], 'k--', lw=0.6)
text(3e-6, 22500, "$n_{\\rm H}~\\widehat{}~{\\tt Cool\\_gamma\\_effective}$", va="top", fontsize=7)
plot([3e-1, 3e1], [15000., 15000.*10.**(2./3.)], 'k--', lw=0.6)
text(3e-1, 200000, "$n_{\\rm H}~\\widehat{}~{\\tt Jeans\\_gamma\\_effective}$", va="top", rotation=43, fontsize=7)
text(0.95e-5, 25, "${\\tt Cool\\_density\\_threshold\\_H\\_p\\_cm3}$", rotation=90, va="bottom", ha="right", fontsize=7)
text(0.95e-1, 25, "${\\tt Jeans\\_density\\_threshold\\_H\\_p\\_cm3}$", rotation=90, va="bottom", ha="right", fontsize=7)
text(5e-8, 8800, "${\\tt Cool\\_temperature\\_norm\\_K}$", va="bottom", fontsize=7)
text(5e-8, 4400, "${\\tt Jeans\\_temperature\\_norm\\_K}$", va="bottom", fontsize=7)
fill_between([1e-5, 1e5], [10, 10], [8000, 8000], color='0.9')
fill_between([1e-1, 1e5], [4000, 400000], color='0.9')
scatter([1e-5], [8000], s=4, color='k')
scatter([1e-1], [4000], s=4, color='k')
xlabel("${\\rm Density}~n_{\\rm H}~[{\\rm cm^{-3}}]$", labelpad=0)
ylabel("${\\rm Temperature}~T~[{\\rm K}]$", labelpad=2)
xlim(3e-8, 3e3)
ylim(20., 2e5)
savefig("EAGLE_entropy_floor.png", dpi=200)
