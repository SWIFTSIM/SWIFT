###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

# Computes the analytical solution of the Toro (2009) test 2 and plots the SPH answer


# Generates the analytical  solution for the Toro (2009) test case 2
# The script works for a given left (x<0) and right (x>0) state and computes the solution at a later time t.
# This follows the solution given in (Toro, 2009)

# Parameters
gas_gamma = 5.0 / 3.0  # Polytropic index
rho_L = 1.0  # Density left state
rho_R = 1.0  # Density right state
v_L = -2.0  # Velocity left state
v_R = 2.0  # Velocity right state
P_L = 0.4  # Pressure left state
P_R = 0.4  # Pressure right state

import sys

sys.path.append("../")
from riemannSolver import RiemannSolver

import matplotlib

matplotlib.use("Agg")
from pylab import *
from scipy import stats
import h5py

style.use("../../../tools/stylesheets/mnras.mplstyle")

snap = int(sys.argv[1])


# Read the simulation data
sim = h5py.File("toroTest2_%04d.hdf5" % snap, "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"]

x = sim["/PartType0/Coordinates"][:, 0]
v = sim["/PartType0/Velocities"][:, 0]
u = sim["/PartType0/InternalEnergies"][:]
S = sim["/PartType0/Entropies"][:]
P = sim["/PartType0/Pressures"][:]
rho = sim["/PartType0/Densities"][:]

x_min = -1.0
x_max = 1.0
x += x_min
N = 1000

# Bin the data
x_bin_edge = np.arange(x_min, x_max, 0.02)
x_bin = 0.5 * (x_bin_edge[1:] + x_bin_edge[:-1])
rho_bin, _, _ = stats.binned_statistic(x, rho, statistic="mean", bins=x_bin_edge)
v_bin, _, _ = stats.binned_statistic(x, v, statistic="mean", bins=x_bin_edge)
P_bin, _, _ = stats.binned_statistic(x, P, statistic="mean", bins=x_bin_edge)
S_bin, _, _ = stats.binned_statistic(x, S, statistic="mean", bins=x_bin_edge)
u_bin, _, _ = stats.binned_statistic(x, u, statistic="mean", bins=x_bin_edge)
rho2_bin, _, _ = stats.binned_statistic(x, rho ** 2, statistic="mean", bins=x_bin_edge)
v2_bin, _, _ = stats.binned_statistic(x, v ** 2, statistic="mean", bins=x_bin_edge)
P2_bin, _, _ = stats.binned_statistic(x, P ** 2, statistic="mean", bins=x_bin_edge)
S2_bin, _, _ = stats.binned_statistic(x, S ** 2, statistic="mean", bins=x_bin_edge)
u2_bin, _, _ = stats.binned_statistic(x, u ** 2, statistic="mean", bins=x_bin_edge)
rho_sigma_bin = np.sqrt(rho2_bin - rho_bin ** 2)
v_sigma_bin = np.sqrt(v2_bin - v_bin ** 2)
P_sigma_bin = np.sqrt(P2_bin - P_bin ** 2)
S_sigma_bin = np.sqrt(S2_bin - S_bin ** 2)
u_sigma_bin = np.sqrt(u2_bin - u_bin ** 2)

# Prepare reference solution
solver = RiemannSolver(gas_gamma)

delta_x = (x_max - x_min) / N
x_s = arange(0.5 * x_min, 0.5 * x_max, delta_x)
rho_s, v_s, P_s, _ = solver.solve(rho_L, v_L, P_L, rho_R, v_R, P_R, x_s / time)
rho_s2, v_s2, P_s2, _ = solver.solve(rho_R, v_R, P_R, rho_L, v_L, P_L, x_s / time)
x_s2 = np.array(x_s)
x_s2 += 1.0
s2neg = x_s2 > 1.0
s2pos = ~s2neg
x_s2[s2neg] -= 2.0

# Additional arrays
u_s = P_s / (rho_s * (gas_gamma - 1.0))  # internal energy
s_s = P_s / rho_s ** gas_gamma  # entropic function
u_s2 = P_s2 / (rho_s2 * (gas_gamma - 1.0))  # internal energy
s_s2 = P_s2 / rho_s2 ** gas_gamma  # entropic function

# Plot the interesting quantities
figure(figsize=(7, 7 / 1.6))

# Velocity profile --------------------------------
subplot(231)
plot(x, v, ".", color="r", ms=0.5, alpha=0.2)
plot(x_s, v_s, "--", color="k", alpha=0.8, lw=1.2)
plot(x_s2[s2pos], v_s2[s2pos], "--", color="k", alpha=0.8, lw=1.2)
plot(x_s2[s2neg], v_s2[s2neg], "--", color="k", alpha=0.8, lw=1.2)
errorbar(x_bin, v_bin, yerr=v_sigma_bin, fmt=".", ms=8.0, color="b", lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Velocity}}~v_x$", labelpad=0)

# Density profile --------------------------------
subplot(232)
plot(x, rho, ".", color="r", ms=0.5, alpha=0.2)
plot(x_s, rho_s, "--", color="k", alpha=0.8, lw=1.2)
plot(x_s2[s2pos], rho_s2[s2pos], "--", color="k", alpha=0.8, lw=1.2)
plot(x_s2[s2neg], rho_s2[s2neg], "--", color="k", alpha=0.8, lw=1.2)
errorbar(x_bin, rho_bin, yerr=rho_sigma_bin, fmt=".", ms=8.0, color="b", lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Density}}~\\rho$", labelpad=0)

# Pressure profile --------------------------------
subplot(233)
plot(x, P, ".", color="r", ms=0.5, alpha=0.2)
plot(x_s, P_s, "--", color="k", alpha=0.8, lw=1.2)
plot(x_s2[s2pos], P_s2[s2pos], "--", color="k", alpha=0.8, lw=1.2)
plot(x_s2[s2neg], P_s2[s2neg], "--", color="k", alpha=0.8, lw=1.2)
errorbar(x_bin, P_bin, yerr=P_sigma_bin, fmt=".", ms=8.0, color="b", lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Pressure}}~P$", labelpad=0)

# Internal energy profile -------------------------
subplot(234)
plot(x, u, ".", color="r", ms=0.5, alpha=0.2)
plot(x_s, u_s, "--", color="k", alpha=0.8, lw=1.2)
plot(x_s2[s2pos], u_s2[s2pos], "--", color="k", alpha=0.8, lw=1.2)
plot(x_s2[s2neg], u_s2[s2neg], "--", color="k", alpha=0.8, lw=1.2)
errorbar(x_bin, u_bin, yerr=u_sigma_bin, fmt=".", ms=8.0, color="b", lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Internal~Energy}}~u$", labelpad=0)

# Entropy profile ---------------------------------
subplot(235)
plot(x, S, ".", color="r", ms=0.5, alpha=0.2)
plot(x_s, s_s, "--", color="k", alpha=0.8, lw=1.2)
plot(x_s2[s2pos], s_s2[s2pos], "--", color="k", alpha=0.8, lw=1.2)
plot(x_s2[s2neg], s_s2[s2neg], "--", color="k", alpha=0.8, lw=1.2)
errorbar(x_bin, S_bin, yerr=S_sigma_bin, fmt=".", ms=8.0, color="b", lw=1.2)
xlabel("${\\rm{Position}}~x$", labelpad=0)
ylabel("${\\rm{Entropy}}~S$", labelpad=0)

# Information -------------------------------------
subplot(236, frameon=False)

text_fontsize = 5

text(
    -0.49,
    0.9,
    "Toro (2009) test 2 with  $\\gamma=%.3f$ in 3D at $t=%.2f$" % (gas_gamma, time),
    fontsize=text_fontsize,
)
text(
    -0.49,
    0.8,
    "Left: $(P_L, \\rho_L, v_L) = (%.3f, %.3f, %.3f)$" % (P_L, rho_L, v_L),
    fontsize=text_fontsize,
)
text(
    -0.49,
    0.7,
    "Right: $(P_R, \\rho_R, v_R) = (%.3f, %.3f, %.3f)$" % (P_R, rho_R, v_R),
    fontsize=text_fontsize,
)
plot([-0.49, 0.1], [0.62, 0.62], "k-", lw=1)
text(-0.49, 0.5, "SWIFT %s" % git.decode("utf-8"), fontsize=text_fontsize)
text(-0.49, 0.4, scheme.decode("utf-8"), fontsize=text_fontsize)
text(-0.49, 0.3, kernel.decode("utf-8"), fontsize=text_fontsize)
text(
    -0.49,
    0.2,
    "$%.2f$ neighbours ($\\eta=%.3f$)" % (neighbours, eta),
    fontsize=text_fontsize,
)
xlim(-0.5, 0.5)
ylim(0, 1)
xticks([])
yticks([])

tight_layout()
savefig("ToroTest2.png", dpi=200)
