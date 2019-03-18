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
    "figure.subplot.left": 0.17,
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

# Metal dependance parameters
Z_0 = 0.002
norm = 0.1
slope = -0.64
max_n = 10.

# Function
Z = logspace(-8, 3, 1000)
n = norm * (Z / Z_0)**slope
n = np.minimum(n, np.ones(np.size(n)) * (max_n))


# Plot the phase space diagram
figure()
subplot(111, xscale="log", yscale="log")

plot(Z, n, 'k-', lw=1.)

plot([3e-4, 3e-2], [1., 1.*100.**(slope)], 'k--', lw=0.6)
plot([1e-10, 1e10], [max_n, max_n], 'k:', lw=0.6)
plot([Z_0, Z_0], [1e-10, norm], 'k:', lw=0.6)
plot([1e-10, Z_0], [norm, norm], 'k:', lw=0.6)
scatter([Z_0], [norm], s=4, color='k')

#arrow(0.014, 0.00025, 0., 0.001, color='k', lw=0.6)
annotate('', xy=(0.014, 1e-3), xytext=(0.014, 3e-4), arrowprops=dict(facecolor='black', shrink=0., width=0.1, headwidth=3., headlength=5.))

text(3e-3, 0.4, "$Z~\\widehat{}~{\\tt threshold\\_slope}$", va="top", rotation=-42, fontsize=7)
text(3e-5, 12., "${\\tt threshold\\_max\\_density\\_H\\_p\\_cm3}$", fontsize=7)
text(3e-7, 0.12, "${\\tt threshold\\_Z0}$", fontsize=7)
text(0.0018, 0.0004, "${\\tt threshold\\_norm\\_H\\_p\\_cm3}$", rotation=90, va="bottom", ha="right", fontsize=7)

xlabel("${\\rm Metallicity~(metal~mass~fraction)}~Z~[-]$", labelpad=2)
ylabel("${\\rm SF~density~threshold}~n_{\\rm H, thresh}~[{\\rm cm^{-3}}]$", labelpad=-1)

xlim(1e-7, 1.)
ylim(0.0002, 50)
savefig("EAGLE_SF_Z_dep.png", dpi=200)
