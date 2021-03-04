###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016  Matthieu Schaller (matthieu.schaller@durham.ac.uk)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
##############################################################################
import matplotlib

matplotlib.use("Agg")
from pylab import *
from scipy import integrate
from scipy import special
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
from matplotlib.font_manager import FontProperties
import numpy
import math

params = {
    "axes.labelsize": 9,
    "axes.titlesize": 10,
    "font.size": 10,
    "legend.fontsize": 10,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "text.usetex": True,
    "figure.figsize": (3.15, 3.15),
    "figure.subplot.left": 0.14,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.1,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0.0,
    "figure.subplot.hspace": 0.0,
    "lines.markersize": 6,
    "lines.linewidth": 3.0,
  #  "text.latex.unicode": True,
}
rcParams.update(params)
rc("font", **{"family": "sans-serif", "sans-serif": ["Times"]})
colors = ["#4477AA", "#CC6677", "#DDCC77", "#117733"]


# Parameters
r_s = 2.0
r_min = 3e-2
r_max = 1.5e2

# Radius
r = logspace(log10(r_min), log10(r_max), 401)
r_rs = r / r_s

k = logspace(log10(r_min / r_s ** 2), log10(r_max / r_s ** 2), 401)
k_rs = k * r_s

# Newtonian solution
phi_newton = 1.0 / r
phit_newton = 1.0 / k ** 2
force_newton = 1.0 / r ** 2


def my_exp(x):
    return (
        1.0
        + x
        + (x ** 2 / 2.0)
        + (x ** 3 / 6.0)
        + (x ** 4 / 24.0)
        + (x ** 5 / 120.0)
        + (x ** 6 / 720.0)
    )
    # return exp(x)


def term(x):  # 1 / (1 + e^x)
    return 1.0 / (1.0 + exp(x))


def my_term(x):  # 1 / (1 + e^x)
    # return 0.5 - 0.25 * x + (x**3 / 48.) - (x**5 / 480)
    return 1.0 / (1.0 + my_exp(x))


def csch(x):  # hyperbolic cosecant
    return 1.0 / sinh(x)


def sigmoid(x):
    return exp(x) * term(x)


def d_sigmoid(x):
    return exp(x) * term(x) ** 2


def my_sigmoid(x):
    # return my_exp(x) / (my_exp(x) + 1.)
    return my_exp(x) * my_term(x)


def my_d_sigmoid(x):
    # return my_exp(x) / ((my_exp(x) + 1)**2)
    return my_exp(x) * my_term(x) ** 2


def swift_corr(x):
    return 2 * sigmoid(4 * x) - 1


def swift_corr2(x):
    return 2 * my_sigmoid(4 * x) - 1


figure()
x = linspace(-4, 4, 100)
plot(x, special.erf(x), "-", color=colors[2])
plot(x, swift_corr(x), "-", color=colors[3])
plot(x, swift_corr2(x), "-.", color=colors[3])
plot(x, x, "-", color=colors[0])
ylim(-1.1, 1.1)
xlim(-4.1, 4.1)
savefig("temp.pdf")


def alpha(x):
    return 1.0 / (1.0 + exp(x))


# Correction in real space
corr_short_gadget2 = special.erf(r / (2.0 * r_s))
corr_short_swift = swift_corr(r / (2.0 * r_s))
corr_short_swift2 = swift_corr2(r / (2.0 * r_s))
eta_short_gadget2 = special.erfc(r / (2.0 * r_s)) + (
    r / (r_s * math.sqrt(math.pi))
) * exp(-r ** 2 / (4.0 * r_s ** 2))
eta_short_swift = (
    4.0 * (r / r_s) * d_sigmoid(2.0 * r / r_s) - 2.0 * sigmoid(2 * r / r_s) + 2.0
)
eta_short_swift2 = (
    4.0 * (r / r_s) * my_d_sigmoid(2.0 * r / r_s) - 2.0 * my_sigmoid(2 * r / r_s) + 2.0
)

# x = 2. * r / r_s
# force_corr = 2. * (1. - exp(x) * (alpha(x) - x * alpha(x)**2))
# force_corr = 2. * (1.- x*exp(x)*alpha(x)**2 - exp(x)*alpha(x))
# force_corr = 2. * (x*alpha(x) - x*alpha(x)**2 -exp(x)*alpha(x) + 1)
# force_corr = abs(2 * (1. - exp(x) * alpha(x) + x * exp(2*x)*alpha(x)**2 - x*exp(x)*alpha(x)))
# force_corr = abs(force_corr)

# Corection in Fourier space
corr_long_gadget2 = exp(-k ** 2 * r_s ** 2)
corr_long_swift = math.pi * k * r_s * csch(0.5 * math.pi * r_s * k) / 2.0

# Shortrange term
phi_short_gadget2 = (1.0 / r) * (1.0 - corr_short_gadget2)
phi_short_swift = (1.0 / r) * (1.0 - corr_short_swift)
phi_short_swift2 = (1.0 / r) * (1.0 - corr_short_swift2)
force_short_gadget2 = (1.0 / r ** 2) * eta_short_gadget2
force_short_swift = (1.0 / r ** 2) * eta_short_swift
force_short_swift2 = (1.0 / r ** 2) * eta_short_swift2

# Long-range term
phi_long_gadget2 = (1.0 / r) * corr_short_gadget2
phi_long_swift = (1.0 / r) * corr_short_swift
phit_long_gadget2 = corr_long_gadget2 / k ** 2
phit_long_swift = corr_long_swift / k ** 2


figure()

# Potential
subplot(311, xscale="log", yscale="log")

plot(r_rs, phi_newton, "--", lw=1.4, label="${\\rm Newtonian}$", color=colors[0])
plot(r_rs, phi_short_gadget2, "-", lw=1.4, label="${\\rm Gadget}$", color=colors[2])
plot(r_rs, phi_short_swift, "-", lw=1.4, label="${\\rm SWIFT}$", color=colors[3])
plot(r_rs, phi_short_swift2, ":", lw=1.4, color=colors[3])
plot([1.0, 1.0], [1e-5, 1e5], "k-.", alpha=0.5, lw=0.5)

xlim(1.1 * r_min / r_s, 0.9 * r_max / r_s)
ylim(1.1 / r_max, 0.9 / r_min)
ylabel("$\\varphi_s(r)$", labelpad=-3)

legend(loc="upper right", frameon=True, handletextpad=0.3, handlelength=1.6, fontsize=8)

# Correction
subplot(312, xscale="log", yscale="log")
plot(r_rs, np.ones(np.size(r)), "--", lw=1.4, color=colors[0])
plot(r_rs, 1.0 - corr_short_gadget2, "-", lw=1.4, color=colors[2])
plot(r_rs, 1.0 - corr_short_swift, "-", lw=1.4, color=colors[3])
plot(r_rs, 1.0 - corr_short_swift2, ":", lw=1.4, color=colors[3])
plot(r_rs, np.ones(np.size(r)) * 0.01, "k-.", alpha=0.5, lw=0.5)
plot([1.0, 1.0], [-1e5, 1e5], "k-.", alpha=0.5, lw=0.5)
plot([-1, -1], [-1, -1], "k-", lw=1.2, label="${\\textrm{Exact}~e^x}$")
plot(
    [-1, -1],
    [-1, -1],
    "k:",
    lw=1.2,
    label="${6^\\textrm{th}~\\textrm{order~series}~e^x}$",
)

yticks([1e-2, 1e-1, 1], ["$0.01$", "$0.1$", "$1$"])
xlim(1.1 * r_min / r_s, 0.9 * r_max / r_s)
ylim(3e-3, 1.5)
# ylabel("$\\chi_s(r)$", labelpad=-3)
ylabel("$\\varphi_s(r) \\times r$", labelpad=-2)

legend(
    loc="center left", frameon=False, handletextpad=0.3, handlelength=1.6, fontsize=7
)

# 1 - Correction
subplot(313, xscale="log", yscale="log")
plot(r_rs, corr_short_gadget2, "-", lw=1.4, color=colors[2])
plot(r_rs, corr_short_swift, "-", lw=1.4, color=colors[3])
plot(r_rs, corr_short_swift2, ":", lw=1.4, color=colors[3])

plot([1.0, 1.0], [1e-5, 1e5], "k-.", alpha=0.5, lw=0.5)
plot(r_rs, np.ones(np.size(r)), "k-.", alpha=0.5, lw=0.5)
plot(r_rs, np.ones(np.size(r)) * 0.01, "k-.", alpha=0.5, lw=0.5)

xlim(1.1 * r_min / r_s, 0.9 * r_max / r_s)
ylim(3e-3, 1.5)
# ylabel("$1 - \\chi_s(r)$", labelpad=-2)
ylabel("$1 - \\varphi_s(r) \\times r$", labelpad=-2)
yticks([1e-2, 1e-1, 1], ["$0.01$", "$0.1$", "$1$"])
xlabel("$r / r_s$", labelpad=1)

savefig("potential_short.pdf")

##################################################################################################


# Force
figure()
subplot(311, xscale="log", yscale="log")

plot(r_rs, force_newton, "--", lw=1.4, label="${\\rm Newtonian}$", color=colors[0])
plot(r_rs, force_short_gadget2, "-", lw=1.4, label="${\\rm Gadget}$", color=colors[2])
plot(r_rs, force_short_swift, "-", lw=1.4, label="${\\rm SWIFT}$", color=colors[3])
# plot(r_rs, (1./r**2) * force_corr, '-', lw=1.2, color='r')
plot(r_rs, force_short_swift2, ":", lw=1.4, color=colors[3])
plot([1.0, 1.0], [1e-5, 1e5], "k-.", alpha=0.5, lw=0.5)

xlim(1.1 * r_min / r_s, 0.9 * r_max / r_s)
ylim(1.1 / r_max ** 2, 0.9 / r_min ** 2)
ylabel("$|\\mathbf{f}_s(r)|$", labelpad=-3)
yticks([1e-4, 1e-2, 1e0, 1e2], ["$10^{-4}$", "$10^{-2}$", "$10^{0}$", "$10^{2}$"])

legend(loc="upper right", frameon=True, handletextpad=0.3, handlelength=1.6, fontsize=8)

# Correction
subplot(312, xscale="log", yscale="log")
plot(r_rs, np.ones(np.size(r)), "--", lw=1.4, color=colors[0])
plot(r_rs, eta_short_gadget2, "-", lw=1.4, color=colors[2])
plot(r_rs, eta_short_swift, "-", lw=1.4, color=colors[3])
plot(r_rs, eta_short_swift2, ":", lw=1.4, color=colors[3])
plot(r_rs, np.ones(np.size(r)) * 0.01, "k-.", alpha=0.5, lw=0.5)
plot([1.0, 1.0], [-1e5, 1e5], "k-.", alpha=0.5, lw=0.5)
plot([-1, -1], [-1, -1], "k-", lw=1.2, label="${\\textrm{Exact}~e^x}$")
plot(
    [-1, -1],
    [-1, -1],
    "k:",
    lw=1.2,
    label="${6^\\textrm{th}~\\textrm{order~series}~e^x}$",
)

yticks([1e-2, 1e-1, 1], ["$0.01$", "$0.1$", "$1$"])
xlim(1.1 * r_min / r_s, 0.9 * r_max / r_s)
ylim(3e-3, 1.5)
ylabel("$|\\mathbf{f}_s(r)|\\times r^2$", labelpad=-2)

legend(
    loc="center left", frameon=False, handletextpad=0.3, handlelength=1.6, fontsize=7
)

# 1 - Correction
subplot(313, xscale="log", yscale="log")
plot(r_rs, 1.0 - eta_short_gadget2, "-", lw=1.4, color=colors[2])
plot(r_rs, 1.0 - eta_short_swift, "-", lw=1.4, color=colors[3])
plot(r_rs, 1.0 - eta_short_swift2, ":", lw=1.4, color=colors[3])

plot([1.0, 1.0], [1e-5, 1e5], "k-.", alpha=0.5, lw=0.5)
plot(r_rs, np.ones(np.size(r)), "k-.", alpha=0.5, lw=0.5)
plot(r_rs, np.ones(np.size(r)) * 0.01, "k-.", alpha=0.5, lw=0.5)

xlim(1.1 * r_min / r_s, 0.9 * r_max / r_s)
ylim(3e-3, 1.5)
ylabel("$1 - |\\mathbf{f}_s(r)|\\times r^2$", labelpad=-3)
yticks([1e-2, 1e-1, 1], ["$0.01$", "$0.1$", "$1$"])
xlabel("$r / r_s$", labelpad=1)

savefig("force_short.pdf")

##################################################################################################

figure()
subplot(311, xscale="log", yscale="log")

# Potential
plot(k_rs, phit_newton, "--", lw=1.4, label="${\\rm Newtonian}$", color=colors[0])
plot(k_rs, phit_long_gadget2, "-", lw=1.4, label="${\\rm Gadget}$", color=colors[2])
plot(k_rs, phit_long_swift, "-", lw=1.4, label="${\\rm SWIFT}$", color=colors[3])
plot([1.0, 1.0], [1e-5, 1e5], "k-.", alpha=0.5, lw=0.5)

legend(loc="lower left", frameon=True, handletextpad=0.3, handlelength=1.6, fontsize=8)

xlim(1.1 * r_min / r_s, 0.9 * r_max / r_s)
ylim(1.1 / r_max ** 2, 0.9 / r_min ** 2)
ylabel("$\\tilde{\\varphi_l}(k)$", labelpad=-3)
yticks([1e-4, 1e-2, 1e0, 1e2], ["$10^{-4}$", "$10^{-2}$", "$10^{0}$", "$10^{2}$"])

subplot(312, xscale="log", yscale="log")

# Potential normalized
plot(
    k_rs,
    phit_newton * k ** 2,
    "--",
    lw=1.4,
    label="${\\rm Newtonian}$",
    color=colors[0],
)
plot(
    k_rs,
    phit_long_gadget2 * k ** 2,
    "-",
    lw=1.4,
    label="${\\rm Gadget}$",
    color=colors[2],
)
plot(
    k_rs, phit_long_swift * k ** 2, "-", lw=1.4, label="${\\rm SWIFT}$", color=colors[3]
)
plot([1.0, 1.0], [1e-5, 1e5], "k-.", alpha=0.5, lw=0.5)
plot(r_rs, np.ones(np.size(r)) * 0.01, "k-.", alpha=0.5, lw=0.5)

xlim(1.1 * r_min / r_s, 0.9 * r_max / r_s)
ylim(3e-3, 1.5)
ylabel("$k^2 \\times \\tilde{\\varphi_l}(k)$", labelpad=-3)
yticks([1e-2, 1e-1, 1], ["$0.01$", "$0.1$", "$1$"])

subplot(313, xscale="log", yscale="log")

plot(
    k_rs,
    1.0 - phit_long_gadget2 * k ** 2,
    "-",
    lw=1.4,
    label="${\\rm Gadget}$",
    color=colors[2],
)
plot(
    k_rs,
    1.0 - phit_long_swift * k ** 2,
    "-",
    lw=1.4,
    label="${\\rm SWIFT}$",
    color=colors[3],
)
plot([1.0, 1.0], [1e-5, 1e5], "k-", alpha=0.5, lw=0.5)
plot(r_rs, np.ones(np.size(r)), "k-.", alpha=0.5, lw=0.5)
plot(r_rs, np.ones(np.size(r)) * 0.01, "k-.", alpha=0.5, lw=0.5)

xlim(1.1 * r_min / r_s, 0.9 * r_max / r_s)
ylim(3e-3, 1.5)
ylabel("$1 - k^2 \\times \\tilde{\\varphi_l}(k)$", labelpad=-3)
yticks([1e-2, 1e-1, 1], ["$0.01$", "$0.1$", "$1$"])

xlabel("$k \\times r_s$", labelpad=1)

savefig("potential_long.pdf")
