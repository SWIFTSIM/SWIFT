import matplotlib

matplotlib.use("Agg")
from pylab import *
from scipy import stats

# Plot parameters
params = {
    "axes.labelsize": 10,
    "axes.titlesize": 10,
    "font.size": 9,
    "font.sans-serif": [
        "Computer Modern",
        "Computer Modern Roman",
        "CMU Serif",
        "cmunrm",
        "DejaVu Sans",
    ],
    "mathtext.fontset": "cm",
    "legend.fontsize": 9,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "text.usetex": False,
    "figure.figsize": (3.15, 3.15),
    "lines.markersize": 6,
    "figure.subplot.left": 0.15,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.13,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0.15,
    "figure.subplot.hspace": 0.12,
    "lines.linewidth": 2.0,
}

rcParams.update(params)

# Metal dependance parameters
Z_0 = 0.002
norm = 0.1
slope = -0.64
max_n = 10.0

# Function
Z = logspace(-8, 3, 1000)
n = norm * (Z / Z_0) ** slope
n = np.minimum(n, np.ones(np.size(n)) * (max_n))


# Plot the phase space diagram
figure()
subplot(111, xscale="log", yscale="log")

plot(Z, n, "k-", lw=1.0)

plot([3e-4, 3e-2], [1.0, 1.0 * 100.0 ** (slope)], "k--", lw=0.6)
plot([1e-10, 1e10], [max_n, max_n], "k:", lw=0.6)
plot([Z_0, Z_0], [1e-10, norm], "k:", lw=0.6)
plot([1e-10, Z_0], [norm, norm], "k:", lw=0.6)
scatter([Z_0], [norm], s=4, color="k")

annotate(
    "",
    xy=(0.014, 1e-3),
    xytext=(0.014, 3e-4),
    arrowprops=dict(
        facecolor="black", shrink=0.0, width=0.1, headwidth=3.0, headlength=5.0
    ),
)
text(0.016, 3.5e-4, "${Z_\\odot}$", fontsize=9)

text(
    3e-4,
    1.45,
    "Z^threshold_slope",
    va="top",
    rotation=-40,
    fontsize=7,
    family="monospace",
)
text(3e-5, 12.0, "threshold_max_density_H_p_cm3", fontsize=7, family="monospace")
text(3e-7, 0.12, "threshold_norm_H_p_cm3", fontsize=7, family="monospace")
text(
    0.0018,
    0.0004,
    "threshold_Z0",
    rotation=90,
    va="bottom",
    ha="right",
    fontsize=7,
    family="monospace",
)

xlabel("Metallicity (metal mass fraction) $Z$ [-]", labelpad=2)
ylabel("SF threshold number density $n_{\\rm H, thresh}$ [cm$^{-3}$]", labelpad=-1)

xlim(1e-7, 1.0)
ylim(0.0002, 50)

savefig("EAGLE_SF_Z_dep.svg", dpi=200)
