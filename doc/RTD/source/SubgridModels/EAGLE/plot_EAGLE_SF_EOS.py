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

# Equations of state
eos_SF_rho = np.logspace(-10, 5, 1000)
eos_SF_T = (eos_SF_rho / 10 ** (-1)) ** (1.0 / 3.0) * 8000.0

# Plot the phase space diagram
figure()
subplot(111, xscale="log", yscale="log")
plot(eos_SF_rho, eos_SF_T, "k-", lw=1.0)

plot([1e-10, 1e-1], [8000, 8000], "k:", lw=0.6)
plot([1e-1, 1e-1], [20, 8000], "k:", lw=0.6)
plot([1e-1, 1e1], [20000.0, 20000.0 * 10.0 ** (2.0 / 3.0)], "k--", lw=0.6)
text(
    0.5e-1,
    200000,
    "$n_{\\rm H}$^EOS_gamma_effective",
    va="top",
    rotation=43,
    fontsize=6.5,
    family="monospace",
)
text(
    0.95e-1,
    25,
    "EOS_density_norm_H_p_cm3",
    rotation=90,
    va="bottom",
    ha="right",
    fontsize=7,
    family="monospace",
)
text(5e-8, 8400, "EOS_temperature_norm_K", va="bottom", fontsize=7, family="monospace")

scatter([1e-1], [8000], s=4, color="k")

xlabel("Hydrogen number density $n_{\\rm H}$ [cm$^{-3}$]", labelpad=0)
ylabel("Temperature $T$ [K]", labelpad=2)
xlim(3e-8, 3e3)
ylim(20.0, 2e5)


savefig("EAGLE_SF_EOS.svg", dpi=200)
