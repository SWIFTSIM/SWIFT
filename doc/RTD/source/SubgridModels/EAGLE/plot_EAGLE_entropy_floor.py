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
eos_cool_rho = np.logspace(-5, 5, 1000)
eos_cool_T = eos_cool_rho ** 0.0 * 8000.0
eos_Jeans_rho = np.logspace(-1, 5, 1000)
eos_Jeans_T = (eos_Jeans_rho / 10 ** (-1)) ** (1.0 / 3.0) * 4000.0

# Plot the phase space diagram
figure()
subplot(111, xscale="log", yscale="log")
plot(eos_cool_rho, eos_cool_T, "k-", lw=1.0)
plot(eos_Jeans_rho, eos_Jeans_T, "k-", lw=1.0)
plot([1e-10, 1e-5], [8000, 8000], "k:", lw=0.6)
plot([1e-10, 1e-1], [4000, 4000], "k:", lw=0.6)
plot([1e-5, 1e-5], [20, 8000], "k:", lw=0.6)
plot([1e-1, 1e-1], [20, 4000], "k:", lw=0.6)
plot([3e-6, 3e-4], [28000, 28000], "k--", lw=0.6)

text(
    3e-6,
    22500,
    "$n_{\\rm H}$^Cool_gamma_effective",
    va="top",
    fontsize=6,
    family="monospace",
)
plot([3e-1, 3e1], [15000.0, 15000.0 * 10.0 ** (2.0 / 3.0)], "k--", lw=0.6)
text(
    3e-1,
    200000,
    "$n_{\\rm H}$^Jeans_gamma_effective",
    va="top",
    rotation=43,
    fontsize=6,
    family="monospace",
)
text(
    0.95e-5,
    23,
    "Cool_density_threshold_H_p_cm3",
    rotation=90,
    va="bottom",
    ha="right",
    fontsize=6,
    family="monospace",
)
text(
    0.95e-1,
    23,
    "Jeans_density_threshold_H_p_cm3",
    rotation=90,
    va="bottom",
    ha="right",
    fontsize=5.5,
    family="monospace",
)
text(5e-8, 8800, "Cool_temperature_norm_K", va="bottom", fontsize=6, family="monospace")
text(
    5e-8, 4400, "Jeans_temperature_norm_K", va="bottom", fontsize=6, family="monospace"
)
fill_between([1e-5, 1e5], [10, 10], [8000, 8000], color="0.9")
fill_between([1e-1, 1e5], [4000, 400000], color="0.9")
scatter([1e-5], [8000], s=4, color="k")
scatter([1e-1], [4000], s=4, color="k")
xlabel("Hydrogen number density $n_{\\rm H}$ [cm$^{-3}$]", labelpad=0)
ylabel("Temperature $T$ [K]", labelpad=2)
xlim(3e-8, 3e3)
ylim(20.0, 2e5)

savefig("EAGLE_entropy_floor.svg", dpi=200)
