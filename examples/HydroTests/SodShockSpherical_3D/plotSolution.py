###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
#               2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
################################################################################

# Compares the swift result for the 2D spherical Sod shock with a high
# resolution 2D reference result

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats
import h5py
import sys

# Parameters
gas_gamma = 5.0 / 3.0  # Polytropic index
rho_L = 1.0  # Density left state
rho_R = 0.125  # Density right state
v_L = 0.0  # Velocity left state
v_R = 0.0  # Velocity right state
P_L = 1.0  # Pressure left state
P_R = 0.1  # Pressure right state

plt.style.use("../../../tools/stylesheets/mnras.mplstyle")

snap = int(sys.argv[1])

# Read the simulation data
sim = h5py.File("sodShock_%04d.hdf5" % snap, "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"]

coords = sim["/PartType0/Coordinates"]
x = np.sqrt(
    (coords[:, 0] - 0.5) ** 2 + (coords[:, 1] - 0.5) ** 2 + (coords[:, 2] - 0.5) ** 2
)
vels = sim["/PartType0/Velocities"]
v = np.sqrt(vels[:, 0] ** 2 + vels[:, 1] ** 2 + vels[:, 2] ** 2)
u = sim["/PartType0/InternalEnergies"][:]
S = sim["/PartType0/Entropies"][:]
P = sim["/PartType0/Pressures"][:]
rho = sim["/PartType0/Densities"][:]

# Bin the data
rho_bin, x_bin_edge, _ = stats.binned_statistic(x, rho, statistic="mean", bins=100)
x_bin = 0.5 * (x_bin_edge[1:] + x_bin_edge[:-1])
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

ref = np.loadtxt("sodShockSpherical3D_exact.txt")

# Plot the interesting quantities
plt.figure(figsize=(7, 7 / 1.6))

line_color = "C4"
binned_color = "C2"
binned_marker_size = 4

scatter_props = dict(
    marker=".",
    ms=1,
    markeredgecolor="none",
    alpha=0.5,
    zorder=-1,
    rasterized=True,
    linestyle="none",
)

errorbar_props = dict(color=binned_color, ms=binned_marker_size, fmt=".", lw=1.2)

# Velocity profile --------------------------------
plt.subplot(231)
plt.plot(x, v, **scatter_props)
plt.plot(ref[:, 0], ref[:, 2], color=line_color, alpha=0.8, lw=1.2)
plt.errorbar(x_bin, v_bin, yerr=v_sigma_bin, **errorbar_props)
plt.xlabel("${\\rm{Radius}}~r$", labelpad=0)
plt.ylabel("${\\rm{Velocity}}~v_r$", labelpad=0)

# Density profile --------------------------------
plt.subplot(232)
plt.plot(x, rho, **scatter_props)
plt.plot(ref[:, 0], ref[:, 1], color=line_color, alpha=0.8, lw=1.2)
plt.errorbar(x_bin, rho_bin, yerr=rho_sigma_bin, **errorbar_props)
plt.xlabel("${\\rm{Radius}}~r$", labelpad=0)
plt.ylabel("${\\rm{Density}}~\\rho$", labelpad=0)

# Pressure profile --------------------------------
plt.subplot(233)
plt.plot(x, P, **scatter_props)
plt.plot(ref[:, 0], ref[:, 3], color=line_color, alpha=0.8, lw=1.2)
plt.errorbar(x_bin, P_bin, yerr=P_sigma_bin, **errorbar_props)
plt.xlabel("${\\rm{Radius}}~r$", labelpad=0)
plt.ylabel("${\\rm{Pressure}}~P$", labelpad=0)

# Internal energy profile -------------------------
plt.subplot(234)
plt.plot(x, u, **scatter_props)
plt.plot(
    ref[:, 0],
    ref[:, 3] / ref[:, 1] / (gas_gamma - 1.0),
    color=line_color,
    alpha=0.8,
    lw=1.2,
)
plt.errorbar(x_bin, u_bin, yerr=u_sigma_bin, **errorbar_props)
plt.xlabel("${\\rm{Radius}}~r$", labelpad=0)
plt.ylabel("${\\rm{Internal~Energy}}~u$", labelpad=0)

# Entropy profile ---------------------------------
plt.subplot(235)
plt.plot(x, S, **scatter_props)
plt.plot(
    ref[:, 0], ref[:, 3] / ref[:, 1] ** gas_gamma, color=line_color, alpha=0.8, lw=1.2
)
plt.errorbar(x_bin, S_bin, yerr=S_sigma_bin, **errorbar_props)
plt.xlabel("${\\rm{Radius}}~r$", labelpad=0)
plt.ylabel("${\\rm{Entropy}}~S$", labelpad=0)

# Information -------------------------------------
plt.subplot(236, frameon=False)

text_fontsize = 5

plt.text(
    -0.45,
    0.9,
    "Sod shock with  $\\gamma=%.3f$ in 3D at $t=%.2f$" % (gas_gamma, time),
    fontsize=text_fontsize,
)
plt.text(
    -0.45,
    0.8,
    "Left:~~ $(P_L, \\rho_L, v_L) = (%.3f, %.3f, %.3f)$" % (P_L, rho_L, v_L),
    fontsize=text_fontsize,
)
plt.text(
    -0.45,
    0.7,
    "Right: $(P_R, \\rho_R, v_R) = (%.3f, %.3f, %.3f)$" % (P_R, rho_R, v_R),
    fontsize=text_fontsize,
)
plt.plot([-0.45, 0.1], [0.62, 0.62], "k-", lw=1)
plt.text(-0.45, 0.5, "$SWIFT$ %s" % git.decode("utf-8"), fontsize=text_fontsize)
plt.text(-0.45, 0.4, scheme.decode("utf-8"), fontsize=text_fontsize)
plt.text(-0.45, 0.3, kernel.decode("utf-8"), fontsize=text_fontsize)
plt.text(
    -0.45,
    0.2,
    "$%.2f$ neighbours ($\\eta=%.3f$)" % (neighbours, eta),
    fontsize=text_fontsize,
)
plt.xlim(-0.5, 0.5)
plt.ylim(0, 1)
plt.xticks([])
plt.yticks([])

plt.tight_layout()

plt.savefig("SodShock.png", dpi=200)
