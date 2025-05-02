###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2023 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
##############################################################################

"""Plot the particle positions from the DemoImpactInitCond settling simulations."""

import os
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import h5py
import woma

# Number of particles
N = 10 ** 5
N_label = "n%d" % (10 * np.log10(N))

# Plotting options
font_size = 20
params = {
    "axes.labelsize": font_size,
    "font.size": font_size,
    "xtick.labelsize": font_size,
    "ytick.labelsize": font_size,
    "font.family": "serif",
}
matplotlib.rcParams.update(params)

# Material colours
Di_mat_colour = {"ANEOS_Fe85Si15": "orangered", "ANEOS_forsterite": "gold"}
Di_id_colour = {woma.Di_mat_id[mat]: colour for mat, colour in Di_mat_colour.items()}

# Scale point size with resolution
size = (1 * np.cbrt(10 ** 6 / N)) ** 2


def load_snapshot(filename):
    """Load and convert the particle data to plot."""
    with h5py.File(filename, "r") as f:
        # Units from file metadata
        file_to_SI = woma.Conversions(
            m=float(f["Units"].attrs["Unit mass in cgs (U_M)"]) * 1e-3,
            l=float(f["Units"].attrs["Unit length in cgs (U_L)"]) * 1e-2,
            t=float(f["Units"].attrs["Unit time in cgs (U_t)"]),
        )

        # Particle data
        A2_pos = (
            np.array(f["PartType0/Coordinates"][()])
            - 0.5 * f["Header"].attrs["BoxSize"]
        ) * file_to_SI.l
        A1_u = np.array(f["PartType0/InternalEnergies"][()]) * file_to_SI.u

    # Restrict to z < 0 for plotting
    A1_sel = np.where(A2_pos[:, 2] < 0)[0]
    A2_pos = A2_pos[A1_sel]
    A1_u = A1_u[A1_sel]

    return A2_pos, A1_u


def plot_snapshot(A2_pos, A1_u):
    """Plot the particles, coloured by their internal energy."""
    plt.figure(figsize=(7, 7))
    ax = plt.gca()
    ax.set_aspect("equal")
    cax = make_axes_locatable(ax).append_axes("right", size="5%", pad=0.05)

    # Earth units
    R_E = 6.3710e6  # m

    # Plot
    scat = ax.scatter(
        A2_pos[:, 0] / R_E,
        A2_pos[:, 1] / R_E,
        c=A1_u,
        edgecolors="none",
        marker=".",
        s=size,
    )
    cbar = plt.colorbar(scat, cax=cax)
    cbar.set_label(r"Sp. Int. Energy (J kg$^{-1}$)")

    ax_lim = 1.5
    ax.set_xlim(-ax_lim, ax_lim)
    ax.set_yticks(ax.get_xticks())
    ax.set_ylim(-ax_lim, ax_lim)
    ax.set_xlabel(r"$x$ ($R_\oplus$)")
    ax.set_ylabel(r"$y$ ($R_\oplus$)")

    plt.tight_layout()


if __name__ == "__main__":
    # Plot each snapshot
    for body in ["target", "impactor"]:
        # Load the data
        snapshot_id = 5
        A2_pos, A1_u = load_snapshot(
            "snapshots/demo_%s_%s_%04d.hdf5" % (body, N_label, snapshot_id)
        )

        # Plot the data
        plot_snapshot(A2_pos, A1_u)

        # Save the figure
        save = "demo_%s_%s_%04d.png" % (body, N_label, snapshot_id)
        plt.savefig(save, dpi=200)
        plt.close()

        print("\rSaved %s" % save)
