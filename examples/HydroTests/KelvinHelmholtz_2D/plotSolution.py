###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016  Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

# Computes the analytical solution of the Gresho-Chan vortex and plots the SPH answer

# Parameters
gas_gamma = 5.0 / 3.0  # Gas adiabatic index
P1 = 2.5  # Central region pressure
P2 = 2.5  # Outskirts pressure
v1 = 0.5  # Central region velocity
v2 = -0.5  # Outskirts vlocity
rho1 = 2  # Central density
rho2 = 1  # Outskirts density

# ---------------------------------------------------------------
# Don't touch anything after this.
# ---------------------------------------------------------------

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import h5py
import sys

plt.style.use("../../../tools/stylesheets/mnras.mplstyle")

snap = int(sys.argv[1])

# Read the simulation data
sim = h5py.File("kelvinHelmholtz_%04d.hdf5" % snap, "r")
boxSize = sim["/Header"].attrs["BoxSize"][0]
time = sim["/Header"].attrs["Time"][0]
scheme = sim["/HydroScheme"].attrs["Scheme"]
kernel = sim["/HydroScheme"].attrs["Kernel function"]
neighbours = sim["/HydroScheme"].attrs["Kernel target N_ngb"]
eta = sim["/HydroScheme"].attrs["Kernel eta"]
git = sim["Code"].attrs["Git Revision"]

pos = sim["/PartType0/Coordinates"][:, :]
x = pos[:, 0] - boxSize / 2
y = pos[:, 1] - boxSize / 2
vel = sim["/PartType0/Velocities"][:, :]
v_norm = np.sqrt(vel[:, 0] ** 2 + vel[:, 1] ** 2)
rho = sim["/PartType0/Densities"][:]
u = sim["/PartType0/InternalEnergies"][:]
S = sim["/PartType0/Entropies"][:]
P = sim["/PartType0/Pressures"][:]

# Plot the interesting quantities
plt.figure(figsize=(7, 7 / 1.6))


# Azimuthal velocity profile -----------------------------
plt.subplot(231)
plt.scatter(
    pos[:, 0],
    pos[:, 1],
    c=vel[:, 0],
    cmap="PuBu",
    edgecolors="face",
    s=0.25,
    vmin=-1.0,
    vmax=1.0,
)
plt.text(
    0.97, 0.97, "${\\rm{Velocity~along}}~x$", ha="right", va="top", backgroundcolor="w"
)
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Position}}~y$", labelpad=0)
plt.xlim(0, 1)
plt.ylim(0, 1)

# Radial density profile --------------------------------
plt.subplot(232)
plt.scatter(
    pos[:, 0],
    pos[:, 1],
    c=rho,
    cmap="PuBu",
    edgecolors="face",
    s=0.25,
    vmin=0.8,
    vmax=2.2,
)
plt.text(0.97, 0.97, "${\\rm{Density}}$", ha="right", va="top", backgroundcolor="w")
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Position}}~y$", labelpad=0)
plt.xlim(0, 1)
plt.ylim(0, 1)

# Radial pressure profile --------------------------------
plt.subplot(233)
plt.scatter(
    pos[:, 0], pos[:, 1], c=P, cmap="PuBu", edgecolors="face", s=0.25, vmin=1, vmax=4
)
plt.text(0.97, 0.97, "${\\rm{Pressure}}$", ha="right", va="top", backgroundcolor="w")
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Position}}~y$", labelpad=0)
plt.xlim(0, 1)
plt.ylim(0, 1)

# Internal energy profile --------------------------------
plt.subplot(234)
plt.scatter(
    pos[:, 0],
    pos[:, 1],
    c=u,
    cmap="PuBu",
    edgecolors="face",
    s=0.25,
    vmin=1.5,
    vmax=5.0,
)
plt.text(
    0.97, 0.97, "${\\rm{Internal~energy}}$", ha="right", va="top", backgroundcolor="w"
)
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Position}}~y$", labelpad=0)
plt.xlim(0, 1)
plt.ylim(0, 1)

# Radial entropy profile --------------------------------
plt.subplot(235)
plt.scatter(
    pos[:, 0],
    pos[:, 1],
    c=S,
    cmap="PuBu",
    edgecolors="face",
    s=0.25,
    vmin=0.5,
    vmax=3.0,
)
plt.text(0.97, 0.97, "${\\rm{Entropy}}$", ha="right", va="top", backgroundcolor="w")
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Position}}~y$", labelpad=0)
plt.xlim(0, 1)
plt.ylim(0, 1)

# Information -------------------------------------
plt.subplot(236, frameon=False)

plt.text(-0.45, 0.9, "Kelvin-Helmholtz instability at $t=%.2f$" % (time), fontsize=10)
plt.text(
    -0.45,
    0.8,
    "Centre: $(P, \\rho, v) = (%.3f, %.3f, %.3f)$" % (P1, rho1, v1),
    fontsize=10,
)
plt.text(
    -0.45,
    0.7,
    "Outskirts: $(P, \\rho, v) = (%.3f, %.3f, %.3f)$" % (P2, rho2, v2),
    fontsize=10,
)
plt.plot([-0.45, 0.1], [0.62, 0.62], "k-", lw=1)
plt.text(-0.45, 0.5, "$SWIFT$ %s" % git.decode("utf-8"), fontsize=10)
plt.text(-0.45, 0.4, scheme.decode("utf-8"), fontsize=10)
plt.text(-0.45, 0.3, kernel.decode("utf-8"), fontsize=10)
plt.text(
    -0.45, 0.2, "$%.2f$ neighbours ($\\eta=%.3f$)" % (neighbours, eta), fontsize=10
)
plt.xlim(-0.5, 0.5)
plt.ylim(0, 1)
plt.xticks([])
plt.yticks([])

plt.tight_layout()

plt.savefig("KelvinHelmholtz.png", dpi=200)
