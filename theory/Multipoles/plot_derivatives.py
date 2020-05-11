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
    "figure.subplot.left": 0.12,
    "figure.subplot.right": 0.99,
    "figure.subplot.bottom": 0.065,
    "figure.subplot.top": 0.99,
    "figure.subplot.wspace": 0.0,
    "figure.subplot.hspace": 0.0,
    "lines.markersize": 6,
    "lines.linewidth": 3.0,
    "text.latex.unicode": True,
}
rcParams.update(params)
rc("font", **{"family": "sans-serif", "sans-serif": ["Times"]})

# Parameters
r_min = 0.0
r_max = 10.0
r_s = 1.7

# Radius
r = linspace(r_min, r_max, 401)
w = 2 * r / r_s

# Powers of alpha ####################################################
alpha = 1.0 / (1.0 + exp(w))
alpha2 = alpha ** 2
alpha3 = alpha ** 3
alpha4 = alpha ** 4
alpha5 = alpha ** 5
alpha6 = alpha ** 6

figure()
plot(w, alpha, label="$\\alpha^1$")
plot(w, alpha2, label="$\\alpha^2$")
plot(w, alpha3, label="$\\alpha^3$")
plot(w, alpha4, label="$\\alpha^4$")
plot(w, alpha5, label="$\\alpha^5$")
plot(w, alpha6, label="$\\alpha^6$")

xlabel("w", labelpad=-4.0)
ylabel("$\\alpha^n(w)$", labelpad=-4.0)

xlim(0, 7.2)
ylim(0.0, 0.52)

legend(loc="upper right")

savefig("alpha_powers.pdf")


# Derivatives of alpha ###############################################
alpha_1 = -alpha + alpha2
alpha_2 = alpha - 3.0 * alpha2 + 2.0 * alpha3
alpha_3 = -alpha + 7.0 * alpha2 - 12.0 * alpha3 + 6.0 * alpha4
alpha_4 = alpha - 15.0 * alpha2 + 50.0 * alpha3 - 60.0 * alpha4 + 24.0 * alpha5
alpha_5 = (
    -alpha
    + 31.0 * alpha2
    - 180.0 * alpha3
    + 390.0 * alpha4
    - 360.0 * alpha5
    + 120.0 * alpha6
)


figure()
plot(w, alpha, label="$\\alpha^{(0)}$")
plot(w, alpha_1, label="$\\alpha^{(1)}$")
plot(w, alpha_2, label="$\\alpha^{(2)}$")
plot(w, alpha_3, label="$\\alpha^{(3)}$")
plot(w, alpha_4, label="$\\alpha^{(4)}$")
plot(w, alpha_5, label="$\\alpha^{(5)}$")

xlabel("w", labelpad=-4.0)
ylabel("$\\alpha^{(n)}(w)$", labelpad=-5.0)

xlim(0, 7.2)
ylim(-0.26, 0.16)

legend(loc="lower right")

savefig("alpha_derivatives.pdf")


# Derivatives of sigma ###############################################
sigma = exp(w) * alpha
sigma_1 = exp(w) * alpha2
sigma_2 = exp(w) * (2 * alpha3 - alpha2)
sigma_3 = exp(w) * (6 * alpha4 - 6 * alpha3 + alpha2)
sigma_4 = exp(w) * (24 * alpha5 - 36 * alpha4 + 14 * alpha3 - alpha2)
sigma_5 = exp(w) * (120 * alpha6 - 240 * alpha5 + 150 * alpha4 - 30 * alpha3 + alpha2)


figure()
plot(w, sigma, label="$\\sigma^{(0)}$")
plot(w, sigma_1, label="$\\sigma^{(1)}$")
plot(w, sigma_2, label="$\\sigma^{(2)}$")
plot(w, sigma_3, label="$\\sigma^{(3)}$")
plot(w, sigma_4, label="$\\sigma^{(4)}$")
plot(w, sigma_5, label="$\\sigma^{(5)}$")

xlabel("w", labelpad=-4.0)
ylabel("$\\sigma^{(n)}(w)$", labelpad=-5.0)

xlim(0, 7.2)
ylim(-0.22, 1.02)

legend(loc="center right")

savefig("sigma_derivatives.pdf")


# Derivatives of chi ###############################################
c1 = 2 / r_s
c2 = (2 / r_s) ** 2
c3 = (2 / r_s) ** 3
c4 = (2 / r_s) ** 4
c5 = (2 / r_s) ** 5

chi = 2 - 2 * exp(w) * alpha
chi_1 = -2 * c1 * exp(w) * alpha2
chi_2 = -2 * c2 * exp(w) * (2 * alpha3 - alpha2)
chi_3 = -2 * c3 * exp(w) * (6 * alpha4 - 6 * alpha3 + alpha2)
chi_4 = -2 * c4 * exp(w) * (24 * alpha5 - 36 * alpha4 + 14 * alpha3 - alpha2)
chi_5 = (
    -2
    * c5
    * exp(w)
    * (120 * alpha6 - 240 * alpha5 + 150 * alpha4 - 30 * alpha3 + alpha2)
)

figure()
plot(r, chi, label="$\\chi^{(0)}$")
plot(r, chi_1, label="$\\chi^{(1)}$")
plot(r, chi_2, label="$\\chi^{(2)}$")
plot(r, chi_3, label="$\\chi^{(3)}$")
plot(r, chi_4, label="$\\chi^{(4)}$")
plot(r, chi_5, label="$\\chi^{(5)}$")

plot([r_s, r_s], [-10, 10], "k--", lw=1)

xlabel("r", labelpad=-4.0)
ylabel("$\\chi^{(n)}(r,r_s)$", labelpad=-5.0)

xlim(0, 7.2)
ylim(-1.52, 1.02)

legend(loc="lower right")

savefig("chi_derivatives.pdf")
