import matplotlib

matplotlib.use("Agg")
from pylab import *
from scipy import special
import numpy as np
import math


e_plummer = 1.0 / 3.0
box_size = 25000
mesh_size = 64
a_smooth = 1.25
r_cut_ratio = 4.5

####################################################################

params = {
    "axes.labelsize": 9,
    "axes.titlesize": 10,
    "font.size": 10,
    "legend.fontsize": 10,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "figure.figsize": (3.15, 3.15),
    "text.latex.unicode": True,
    "text.usetex": True,
}
rcParams.update(params)

plummer_to_spline_ratio = 3.0

H = plummer_to_spline_ratio * e_plummer
r_s = a_smooth * box_size / mesh_size
r_cut = r_s * r_cut_ratio

MAC_lo_limit = (5.0 / 9.0) * H
MAC_hi_limit = (5.0 / 3.0) * r_s

print(("Potential softened below", H, "kpc and truncated above", r_s, "kpc"))

####################################################################

r = np.logspace(np.log10(e_plummer) - 1.2, np.log10(box_size) + 0.2, 10000)

# Newtonian gravity
f_newton = 1 / r ** 2

# Simulated gravity
u = r / H
u = u[u <= 1]

W_swift = 21.0 * u ** 6 - 90.0 * u ** 5 + 140.0 * u ** 4 - 84.0 * u ** 3 + 14.0 * u
f_swift = f_newton * (
    special.erfc(0.5 * r / r_s)
    + (1.0 / math.sqrt(math.pi)) * (r / r_s) * np.exp(-0.25 * (r / r_s) ** 2)
)
f_swift[r <= H] = W_swift / H ** 2
f_swift[r > r_cut] = 0

W_gadget = u * (
    21.333333 - 48 * u + 38.4 * u ** 2 - 10.6666667 * u ** 3 - 0.06666667 * u ** -3
)
W_gadget[u < 0.5] = u[u < 0.5] * (
    10.666667 + u[u < 0.5] ** 2 * (32.0 * u[u < 0.5] - 38.4)
)
f_gadget = f_newton * (
    special.erfc(0.5 * r / r_s)
    + (1.0 / math.sqrt(math.pi)) * (r / r_s) * np.exp(-0.25 * (r / r_s) ** 2)
)
f_gadget[r <= H] = W_gadget / H ** 2
f_gadget[r > r_cut] = 0

f_MAC = np.copy(f_newton)
f_MAC[r < MAC_lo_limit] = (1 / r[r < MAC_lo_limit]) ** 0 * MAC_lo_limit ** -2
f_MAC[r > MAC_hi_limit] = (1 / r[r > MAC_hi_limit]) ** 4 * MAC_hi_limit ** 2
f_MAC[r > r_cut] = 0

# range_test = np.logical_and(r > 0.01 * e_plummer, r < 2 * r_cut)
# print(np.max(f_swift[range_test] / f_MAC[range_test]))

####################################################################

fig = figure()
colors = ["#4477AA", "#CC6677", "#DDCC77", "#117733"]
gs1 = fig.add_gridspec(
    nrows=4,
    ncols=1,
    left=0.14,
    right=0.99,
    wspace=0.0,
    hspace=0.0,
    top=0.99,
    bottom=0.1,
)
fig.add_subplot(gs1[0:3, :], xscale="log", yscale="log")

plot(r, f_newton, "--", color=colors[0], label="Newtonian")
plot(r, f_swift, "-", color=colors[3], label="SWIFT")
plot(r, f_MAC, "-.", color=colors[2], label="MAC estimator")
# plot(r, f_gadget, '-.', color=colors[2], label="Gadget")

plot([e_plummer, e_plummer], [1e-20, 1e20], "k--", alpha=0.3, lw=0.7)
plot([H, H], [1e-20, 1e20], "k--", alpha=0.3, lw=0.7)
plot([r_s, r_s], [1e-20, 1e20], "k--", alpha=0.3, lw=0.7)
plot([r_cut, r_cut], [1e-20, 1e20], "k--", alpha=0.3, lw=0.7)
plot([box_size, box_size], [1e-20, 1e20], "k--", alpha=0.3, lw=0.7)

text(
    e_plummer,
    1e-9,
    "$\\epsilon_{\\rm Plummer}$",
    rotation=90,
    backgroundcolor="w",
    ha="center",
    alpha=0.3,
)
# text(H, 1e-9, "$\\epsilon_{\\rm spline}$", rotation=90, backgroundcolor='w', ha="center", alpha=0.3)
text(H, 1e-9, "$H$", rotation=90, backgroundcolor="w", ha="center", alpha=0.3)
text(
    r_s,
    1e-1,
    "$r_{\\rm s}$",
    rotation=90,
    backgroundcolor="w",
    ha="center",
    va="top",
    alpha=0.3,
)
text(
    r_cut,
    1e-1,
    "$r_{\\rm cut}$",
    rotation=90,
    backgroundcolor="w",
    ha="center",
    va="top",
    alpha=0.3,
)
text(
    box_size,
    1e-1,
    "$L$",
    rotation=90,
    backgroundcolor="w",
    ha="center",
    va="top",
    alpha=0.3,
)

legend(
    loc="upper right",
    frameon=True,
    handletextpad=0.3,
    handlelength=1.6,
    fontsize=8,
    framealpha=1.0,
)

ylim(0.1 * (box_size) ** -2, 2 * (e_plummer / 30) ** -2)
xlim(e_plummer / 30, box_size * 2.5)

tick_params(axis="x", which="both", labelbottom=False)

xlabel("$r$")
ylabel("$|f(r)|$", labelpad=-2)

##################################################################################
fig.add_subplot(gs1[3, :], xscale="log", yscale="log")


plot(r, f_newton * r ** 2, "--", color=colors[0], label="Newtonian")
plot(r, f_swift * r ** 2, "-", color=colors[3], label="SWIFT")
plot(r, f_MAC * r ** 2, "-.", color=colors[2], label="MAC estimator")
# plot(r, f_gadget * r**2, '-.', color=colors[2], label="Gadget")

plot([e_plummer, e_plummer], [1e-20, 1e20], "k--", alpha=0.3, lw=0.7)
plot([H, H], [1e-20, 1e20], "k--", alpha=0.3, lw=0.7)
plot([r_s, r_s], [1e-20, 1e20], "k--", alpha=0.3, lw=0.7)
plot([r_cut, r_cut], [1e-20, 1e20], "k--", alpha=0.3, lw=0.7)
plot([box_size, box_size], [1e-20, 1e20], "k--", alpha=0.3, lw=0.7)

ylim(0.08, 2.2)
xlim(e_plummer / 30, box_size * 2.5)

yticks([0.1, 1], ["$0.1$", "$1$"])

xlabel("$r$", labelpad=0)
ylabel("$|f(r)| \\times r^2$", labelpad=0)

##################################################################################
# fig.add_subplot(gs1[4, :], xscale="log", yscale="log")

# plot(r, f_newton / f_swift, '--', color=colors[0], label="Newtonian")
# plot(r, f_swift / f_swift, '-', color=colors[3], label="SWIFT")
# plot(r, f_MAC / f_swift, ':', color=colors[1], label="MAC estimator")
# plot(r, f_gadget / f_swift, '-.', color=colors[2], label="Gadget")

# plot([e_plummer, e_plummer], [1e-20, 1e20], 'k--', alpha=0.3, lw=0.7)
# plot([H, H], [1e-20, 1e20], 'k--', alpha=0.3, lw=0.7)
# plot([r_s, r_s], [1e-20, 1e20], 'k--', alpha=0.3, lw=0.7)
# plot([r_cut, r_cut], [1e-20, 1e20], 'k--', alpha=0.3, lw=0.7)
# plot([box_size, box_size], [1e-20, 1e20], 'k--', alpha=0.3, lw=0.7)

# ylim(0.5, 13.)
# xlim(e_plummer / 30, box_size * 1.6)

# xlabel("$r$", labelpad=0)
# ylabel("$|f(r)| / |f_{SWIFT}(r)|$", labelpad=2)


savefig("mac_potential.pdf")
