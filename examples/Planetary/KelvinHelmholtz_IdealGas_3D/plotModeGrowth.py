###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2025 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
#               2025 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
"""
Generate plot Mode growth of the 3D Kelvin--Helmholtz instability.
This is based on the quantity calculated in Eqns. 10--13 of McNally et al. 2012
"""

import sys
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


def calculate_mode_growth(snap, vy_init_amp, wavelength):

    # Load data
    snap_file = "kelvin_helmholtz_%04d.hdf5" % snap

    with h5py.File(snap_file, "r") as f:
        boxsize_l = f["Header"].attrs["BoxSize"][0]
        A2_pos = f["/PartType0/Coordinates"][:, :] / boxsize_l
        A2_vel = f["/PartType0/Velocities"][:, :] / (vy_init_amp * boxsize_l)
        A1_h = f["/PartType0/SmoothingLengths"][:] / boxsize_l

    # Adjust positions
    A2_pos[:, 0] += 0.5
    A2_pos[:, 1] += 0.5

    # Masks to select the upper and lower halfs of the simulations
    mask_up = A2_pos[:, 1] >= 0.5
    mask_down = A2_pos[:, 1] < 0.5

    # McNally et al. 2012 Eqn. 10
    s = np.empty(len(A1_h))
    s[mask_down] = (
        A2_vel[mask_down, 1]
        * A1_h[mask_down] ** 3
        * np.sin(2 * np.pi * A2_pos[mask_down, 0] / wavelength)
        * np.exp(-2 * np.pi * np.abs(A2_pos[mask_down, 1] - 0.25) / wavelength)
    )
    s[mask_up] = (
        A2_vel[mask_up, 1]
        * A1_h[mask_up] ** 3
        * np.sin(2 * np.pi * A2_pos[mask_up, 0] / wavelength)
        * np.exp(-2 * np.pi * np.abs((1 - A2_pos[mask_up, 1]) - 0.25) / wavelength)
    )

    # McNally et al. 2012 Eqn. 11
    c = np.empty(len(A1_h))
    c[mask_down] = (
        A2_vel[mask_down, 1]
        * A1_h[mask_down] ** 3
        * np.cos(2 * np.pi * A2_pos[mask_down, 0] / wavelength)
        * np.exp(-2 * np.pi * np.abs(A2_pos[mask_down, 1] - 0.25) / wavelength)
    )
    c[mask_up] = (
        A2_vel[mask_up, 1]
        * A1_h[mask_up] ** 3
        * np.cos(2 * np.pi * A2_pos[mask_up, 0] / wavelength)
        * np.exp(-2 * np.pi * np.abs((1 - A2_pos[mask_up, 1]) - 0.25) / wavelength)
    )

    # McNally et al. 2012 Eqn. 12
    d = np.empty(len(A1_h))
    d[mask_down] = A1_h[mask_down] ** 3 * np.exp(
        -2 * np.pi * np.abs(A2_pos[mask_down, 1] - 0.25) / wavelength
    )
    d[mask_up] = A1_h[mask_up] ** 3 * np.exp(
        -2 * np.pi * np.abs((1 - A2_pos[mask_up, 1]) - 0.25) / wavelength
    )

    # McNally et al. 2012 Eqn. 13
    M = 2 * np.sqrt((np.sum(s) / np.sum(d)) ** 2 + (np.sum(c) / np.sum(d)) ** 2)

    return M


if __name__ == "__main__":

    # Simulation paramerters for nomralisation of mode growth
    vy_init_amp = (
        0.01
    )  # Initial amplitude of y velocity perturbation in units of the boxsize per second
    wavelength = 0.5  # wavelength of initial perturbation in units of the boxsize

    # Use Gridspec to set up figure
    n_gs_ax = 40
    ax_len = 5
    fig = plt.figure(figsize=(9, 6))
    gs = mpl.gridspec.GridSpec(nrows=40, ncols=60)
    ax = plt.subplot(gs[:, :])

    # Snapshots and corresponding times
    snaps = np.arange(21)
    times = 0.1 * snaps

    # Calculate mode mode growth
    mode_growth = np.empty(len(snaps))
    for i, snap in enumerate(snaps):
        M = calculate_mode_growth(snap, vy_init_amp, wavelength)
        mode_growth[i] = M

    # Plot
    ax.plot(times, np.log10(mode_growth), linewidth=1.5)

    ax.set_ylabel(r"$\log( \, M \; / \; M^{}_0 \, )$", fontsize=18)
    ax.set_xlabel(r"$t \; / \; \tau^{}_{\rm KH}$", fontsize=18)
    ax.minorticks_on()
    ax.tick_params(which="major", direction="in")
    ax.tick_params(which="minor", direction="in")
    ax.tick_params(labelsize=14)

    plt.savefig("mode_growth.png", dpi=300, bbox_inches="tight")
