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
v = sim["/PartType0/Velocities"][:, 0]
u = sim["/PartType0/InternalEnergies"][:]
S = sim["/PartType0/Entropies"][:]
P = sim["/PartType0/Pressures"][:]
rho = sim["/PartType0/Densities"][:]

N = 1001  # Number of points
x_min = -1
x_max = 1

x += x_min

# ---------------------------------------------------------------
# Don't touch anything after this.
# ---------------------------------------------------------------
x_s = np.arange(0, 2.0, 2.0 / N) - 1.0
rho_s = np.ones(N) * rho0
P_s = np.ones(N) * rho0
v_s = np.ones(N) * v0

# Shock position
u0 = rho0 * P0 * (gas_gamma - 1)
us = 0.5 * (gas_gamma - 1) * v0
rs = us * time

# Post-shock values
rho_s[np.abs(x_s) < rs] = rho0 * (gas_gamma + 1) / (gas_gamma - 1)
P_s[np.abs(x_s) < rs] = 0.5 * rho0 * v0 ** 2 * (gas_gamma + 1)
v_s[np.abs(x_s) < rs] = 0.0

# Pre-shock values
rho_s[np.abs(x_s) >= rs] = rho0
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
    ms=4,
    markeredgecolor="none",
    zorder=-1,
    rasterized=True,
    linestyle="none",
)

errorbar_props = dict(color=binned_color, ms=binned_marker_size, fmt=".", lw=1.2)

# Velocity profile --------------------------------
plt.subplot(231)
plt.plot(x, v, **scatter_props)
plt.plot(x_s, v_s, "--", color=line_color, alpha=0.8, lw=1.2)
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Velocity}}~v_x$", labelpad=-4)
plt.xlim(-0.5, 0.5)
plt.ylim(-1.2, 1.2)

# Density profile --------------------------------
plt.subplot(232)
plt.plot(x, rho, **scatter_props)
plt.plot(x_s, rho_s, "--", color=line_color, alpha=0.8, lw=1.2)
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Density}}~\\rho$", labelpad=0)
plt.xlim(-0.5, 0.5)
plt.ylim(0.95, 4.4)

# Pressure profile --------------------------------
plt.subplot(233)
plt.plot(x, P, **scatter_props)
plt.plot(x_s, P_s, "--", color=line_color, alpha=0.8, lw=1.2)
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Pressure}}~P$", labelpad=0)
plt.xlim(-0.5, 0.5)
plt.ylim(-0.05, 1.8)

# Internal energy profile -------------------------
plt.subplot(234)
plt.plot(x, u, **scatter_props)
plt.plot(x_s, u_s, "--", color=line_color, alpha=0.8, lw=1.2)
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Internal~Energy}}~u$", labelpad=0)
plt.xlim(-0.5, 0.5)
plt.ylim(-0.05, 0.8)

# Entropy profile ---------------------------------
plt.subplot(235)
plt.plot(x, S, **scatter_props)
plt.plot(x_s, s_s, "--", color=line_color, alpha=0.8, lw=1.2)
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Entropy}}~S$", labelpad=-9)
plt.xlim(-0.5, 0.5)
plt.ylim(-0.05, 0.2)

# Information -------------------------------------
plt.subplot(236, frameon=False)

text_fontsize = 5

plt.text(
    -0.45,
    0.9,
    "Noh problem with  $\\gamma=%.3f$ in 1D at $t=%.2f$" % (gas_gamma, time),
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

plt.savefig("Noh.png", dpi=200)
