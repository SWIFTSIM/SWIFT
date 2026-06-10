###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

# Computes the analytical solution of the Noh problem and plots the SPH answer


# Parameters
gas_gamma = 5.0 / 3.0  # Polytropic index
rho0 = 1.0  # Background density
P0 = 1.0e-6  # Background pressure
v0 = 1

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy import stats
import numpy as np
import h5py
import sys

plt.style.use("../../../tools/stylesheets/mnras.mplstyle")

snap = int(sys.argv[1])

# Read the simulation data
sim = h5py.File("noh_%04d.hdf5" % snap, "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"]

x = sim["/PartType0/Coordinates"][:, 0]
y = sim["/PartType0/Coordinates"][:, 1]
z = sim["/PartType0/Coordinates"][:, 2]
vx = sim["/PartType0/Velocities"][:, 0]
vy = sim["/PartType0/Velocities"][:, 1]
vz = sim["/PartType0/Velocities"][:, 2]
Bx = sim["/PartType0/MagneticFluxDensities"][:, 0]
By = sim["/PartType0/MagneticFluxDensities"][:, 1]
Bz = sim["/PartType0/MagneticFluxDensities"][:, 2]
u = sim["/PartType0/InternalEnergies"][:]
S = sim["/PartType0/Entropies"][:]
P = sim["/PartType0/Pressures"][:]
rho = sim["/PartType0/Densities"][:]

r = np.sqrt((x - 1) ** 2 + (y - 1) ** 2 + (z - 1) ** 2)
v = -np.sqrt(vx ** 2 + vy ** 2 + vz ** 2)
Pmag = 0.5 * (Bx ** 2 + By ** 2 + Bz ** 2)

# Bin the data
r_bin_edge = np.arange(0.0, 1.0, 0.02)
r_bin = 0.5 * (r_bin_edge[1:] + r_bin_edge[:-1])
rho_bin, _, _ = stats.binned_statistic(r, rho, statistic="mean", bins=r_bin_edge)
v_bin, _, _ = stats.binned_statistic(r, v, statistic="mean", bins=r_bin_edge)
P_bin, _, _ = stats.binned_statistic(r, P, statistic="mean", bins=r_bin_edge)
Pmag_bin, _, _ = stats.binned_statistic(r, Pmag, statistic="mean", bins=r_bin_edge)
S_bin, _, _ = stats.binned_statistic(r, S, statistic="mean", bins=r_bin_edge)
u_bin, _, _ = stats.binned_statistic(r, u, statistic="mean", bins=r_bin_edge)
rho2_bin, _, _ = stats.binned_statistic(r, rho ** 2, statistic="mean", bins=r_bin_edge)
v2_bin, _, _ = stats.binned_statistic(r, v ** 2, statistic="mean", bins=r_bin_edge)
P2_bin, _, _ = stats.binned_statistic(r, P ** 2, statistic="mean", bins=r_bin_edge)
Pmag2_bin, _, _ = stats.binned_statistic(
    r, Pmag ** 2, statistic="mean", bins=r_bin_edge
)
S2_bin, _, _ = stats.binned_statistic(r, S ** 2, statistic="mean", bins=r_bin_edge)
u2_bin, _, _ = stats.binned_statistic(r, u ** 2, statistic="mean", bins=r_bin_edge)
rho_sigma_bin = np.sqrt(rho2_bin - rho_bin ** 2)
v_sigma_bin = np.sqrt(v2_bin - v_bin ** 2)
P_sigma_bin = np.sqrt(P2_bin - P_bin ** 2)
Pmag_sigma_bin = np.sqrt(Pmag2_bin - Pmag_bin ** 2)
S_sigma_bin = np.sqrt(S2_bin - S_bin ** 2)
u_sigma_bin = np.sqrt(u2_bin - u_bin ** 2)


# Analytic solution
N = 1000  # Number of points

x_s = np.arange(0, 2.0, 2.0 / N) - 1.0
rho_s = np.ones(N) * rho0
P_s = np.ones(N) * rho0
v_s = np.ones(N) * v0

# Shock position
u0 = rho0 * P0 * (gas_gamma - 1)
us = 0.5 * (gas_gamma - 1) * v0
rs = us * time

# Post-shock values
rho_s[np.abs(x_s) < rs] = rho0 * ((gas_gamma + 1) / (gas_gamma - 1)) ** 3
P_s[np.abs(x_s) < rs] = (
    0.5 * rho0 * v0 ** 2 * (gas_gamma + 1) ** 3 / (gas_gamma - 1) ** 2
)
v_s[np.abs(x_s) < rs] = 0.0

# Pre-shock values
rho_s[np.abs(x_s) >= rs] = rho0 * (1 + v0 * time / np.abs(x_s[np.abs(x_s) >= rs])) ** 2
P_s[np.abs(x_s) >= rs] = 0
v_s[x_s >= rs] = -v0
v_s[x_s <= -rs] = v0

# Additional arrays
u_s = P_s / (rho_s * (gas_gamma - 1.0))  # internal energy
s_s = P_s / rho_s ** gas_gamma  # entropic function

# Plot the interesting quantities
plt.figure(figsize=(7, 7 / 1.6))

line_color = "C4"
binned_color = "C2"
binned_marker_size = 4

scatter_props = dict(
    marker=".",
    ms=1,
    markeredgecolor="none",
    alpha=0.2,
    zorder=-1,
    rasterized=True,
    linestyle="none",
)

errorbar_props = dict(color=binned_color, ms=binned_marker_size, fmt=".", lw=1.2)

# Velocity profile --------------------------------
plt.subplot(231)
plt.plot(r, v, **scatter_props)
plt.plot(x_s, v_s, "--", color=line_color, alpha=0.8, lw=1.2)
plt.errorbar(r_bin, v_bin, yerr=v_sigma_bin, **errorbar_props)
plt.xlabel("${\\rm{Radius}}~r$", labelpad=0)
plt.ylabel("${\\rm{Velocity}}~v_r$", labelpad=-4)
plt.xlim(0, 0.5)
plt.ylim(-1.2, 0.4)

# Density profile --------------------------------
plt.subplot(232)
plt.plot(r, rho, **scatter_props)
plt.plot(x_s, rho_s, "--", color=line_color, alpha=0.8, lw=1.2)
plt.errorbar(r_bin, rho_bin, yerr=rho_sigma_bin, **errorbar_props)
plt.xlabel("${\\rm{Radius}}~r$", labelpad=0)
plt.ylabel("${\\rm{Density}}~\\rho$", labelpad=0)
plt.xlim(0, 0.5)
plt.ylim(0.95, 71)

# Pressure profile --------------------------------
plt.subplot(233)
plt.plot(r, P, **scatter_props)
plt.plot(x_s, P_s, "--", color=line_color, alpha=0.8, lw=1.2)
plt.errorbar(r_bin, P_bin, yerr=P_sigma_bin, **errorbar_props)
plt.xlabel("${\\rm{Radius}}~r$", labelpad=0)
plt.ylabel("${\\rm{Pressure}}~P$", labelpad=0)
plt.xlim(0, 0.5)
plt.ylim(-0.5, 25)

# Internal energy profile -------------------------
plt.subplot(234)
plt.plot(r, u, **scatter_props)
plt.plot(x_s, u_s, "--", color=line_color, alpha=0.8, lw=1.2)
plt.errorbar(r_bin, u_bin, yerr=u_sigma_bin, **errorbar_props)
plt.xlabel("${\\rm{Radius}}~r$", labelpad=0)
plt.ylabel("${\\rm{Internal~Energy}}~u$", labelpad=0)
plt.xlim(0, 0.5)
plt.ylim(-0.05, 0.8)

"""
# Entropy profile ---------------------------------
plt.subplot(235)
plt.plot(r, S, **scatter_props)
plt.plot(x_s, s_s, "--", color=line_color, alpha=0.8, lw=1.2)
plt.errorbar(r_bin, S_bin, yerr=S_sigma_bin, **errorbar_props)
plt.xlabel("${\\rm{Radius}}~r$", labelpad=0)
plt.ylabel("${\\rm{Entropy}}~S$", labelpad=-9)
plt.xlim(0, 0.5)
plt.ylim(-0.05, 0.2)
"""

# Magnetic Pressure -------------------------------
plt.subplot(235)
plt.semilogy(r, Pmag / (1.256637e01 * u), **scatter_props)
plt.semilogy(r, Pmag / (1.256637e01 * 0.5 * v ** 2), **scatter_props)
# plt.plot(x_s, s_s, "--", color=line_color, alpha=0.8, lw=1.2)
# plt.errorbar(r_bin, Pmag_bin, yerr=Pmag_sigma_bin, **errorbar_props)
plt.xlabel("${\\rm{Radius}}~r$", labelpad=0)
plt.ylabel("${\\rm{Energy Ratios}}~l$", labelpad=-9)
plt.xlim(0, 0.5)
# plt.ylim(-0.05, 0.2)

# Information -------------------------------------
plt.subplot(236, frameon=False)

text_fontsize = 5

plt.text(
    -0.45,
    0.9,
    "Noh problem with  $\\gamma=%.3f$ in 3D at $t=%.2f$" % (gas_gamma, time),
    fontsize=text_fontsize,
)
plt.text(
    -0.45,
    0.8,
    "ICs: $(P_0, \\rho_0, v_0) = (%1.2e, %.3f, %.3f)$" % (1e-6, 1.0, -1.0),
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

plt.savefig(sys.argv[2], dpi=200)
