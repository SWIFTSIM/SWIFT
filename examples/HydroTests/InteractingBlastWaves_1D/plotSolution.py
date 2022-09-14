###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
import matplotlib.pyplot as plt
import numpy as np
import sys
import h5py

plt.style.use("../../../tools/stylesheets/mnras.mplstyle")

snap = int(sys.argv[1])

# Open the file and read the relevant data
file = h5py.File("interactingBlastWaves_{0:04d}.hdf5".format(snap), "r")
x = file["/PartType0/Coordinates"][:, 0]
rho = file["/PartType0/Densities"][:]
v = file["/PartType0/Velocities"][:, 0]
u = file["/PartType0/InternalEnergies"][:]
S = file["/PartType0/Entropies"][:]
P = file["/PartType0/Pressures"][:]
time = file["/Header"].attrs["Time"][0]

scheme = file["/HydroScheme"].attrs["Scheme"]
kernel = file["/HydroScheme"].attrs["Kernel function"]
neighbours = file["/HydroScheme"].attrs["Kernel target N_ngb"][0]
eta = file["/HydroScheme"].attrs["Kernel eta"][0]
gamma = file["/HydroScheme"].attrs["Adiabatic index"]
git = file["Code"].attrs["Git Revision"]

if gamma != 1.4:
    print(
        "Error: SWIFT was run with the wrong adiabatic index. Should have been 1.4",
        gamma,
    )
    exit(1)

ref = np.loadtxt("interactingBlastWaves1D_exact.txt")

# Plot the interesting quantities
plt.figure(figsize=(7, 7 / 1.6))

line_color = "C4"
binned_color = "C2"
binned_marker_size = 4

scatter_props = dict(
    marker=".",
    ms=4,
    markeredgecolor="none",
    alpha=0.2,
    zorder=-1,
    rasterized=True,
    linestyle="none",
)

# Velocity profile
plt.subplot(231)
plt.plot(x, v, **scatter_props)
plt.plot(ref[:, 0], ref[:, 2], "--", color=line_color, alpha=0.8, lw=1.2)
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Velocity}}~v_x$", labelpad=0)
plt.xlim(0.0, 1.0)
plt.ylim(-1.0, 15.0)

# Density profile
plt.subplot(232)
plt.plot(x, rho, **scatter_props)
plt.plot(ref[:, 0], ref[:, 1], "--", color=line_color, alpha=0.8, lw=1.2)
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Density}}~\\rho$", labelpad=0)
plt.xlim(0.0, 1.0)

# Pressure profile
plt.subplot(233)
plt.plot(x, P, **scatter_props)
plt.plot(ref[:, 0], ref[:, 3], "--", color=line_color, alpha=0.8, lw=1.2)
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Pressure}}~P$", labelpad=0)
plt.xlim(0.0, 1.0)

# Internal energy profile
plt.subplot(234)
plt.plot(x, u, **scatter_props)
plt.plot(
    ref[:, 0],
    ref[:, 3] / ref[:, 1] / (gamma - 1.0),
    "--",
    color=line_color,
    alpha=0.8,
    lw=1.2,
)
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Internal~Energy}}~u$", labelpad=0)
plt.xlim(0.0, 1.0)

# Entropy profile
plt.subplot(235)
plt.plot(x, S, **scatter_props)
plt.plot(
    ref[:, 0], ref[:, 3] / ref[:, 1] ** gamma, "--", color=line_color, alpha=0.8, lw=1.2
)
plt.xlabel("${\\rm{Position}}~x$", labelpad=0)
plt.ylabel("${\\rm{Entropy}}~S$", labelpad=0)
plt.xlim(0.0, 1.0)

# Run information
plt.subplot(236, frameon=False)

text_fontsize = 5

plt.text(
    -0.45,
    0.9,
    "Interacting blast waves test\nwith $\\gamma={0:.3f}$ in 1D at $t = {1:.2f}$".format(
        gamma[0], time
    ),
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
plt.ylim(0.0, 1.0)
plt.xticks([])
plt.yticks([])

plt.tight_layout()

plt.savefig("InteractingBlastWaves.png", dpi=200)
