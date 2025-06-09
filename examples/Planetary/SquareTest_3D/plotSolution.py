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
"""Plot a central slice from a 3D square test snapshot.

Parameters
----------
type : str
    Either "equal_spacing" or "equal_mass".

snap : int
    The snapshot ID to plot.
"""

import sys
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


def plot_square(A1_x, A1_y, A1_rho, A1_u, A1_P, A1_size, boxsize_l):
    # Use Gridspec to set up figure
    n_gs_ax = 40
    n_gs_ax_gap = 10
    n_gs_cbar_gap = 2
    n_gs_cbar = 2

    ax_len = 5
    ax_gap_len = n_gs_ax_gap * ax_len / n_gs_ax
    cbar_gap_len = n_gs_cbar_gap * ax_len / n_gs_ax
    cbar_len = n_gs_cbar * ax_len / n_gs_ax

    fig = plt.figure(
        figsize=(3 * ax_len + 2 * ax_gap_len + 3 * cbar_gap_len + 3 * cbar_len, ax_len)
    )
    gs = mpl.gridspec.GridSpec(
        nrows=n_gs_ax,
        ncols=3 * n_gs_ax + 2 * n_gs_ax_gap + 3 * n_gs_cbar_gap + 3 * n_gs_cbar,
    )

    # Quantities to plot
    plot_quantity = ["Density", "Internal Energy", "Pressure"]
    plot_vectors = [A1_rho, A1_u, A1_P]
    cmaps = ["Spectral_r", "inferno", "viridis"]

    for i in range(3):
        # Set up subfig and color bar axes
        y0 = 0
        y1 = n_gs_ax
        x0 = i * (n_gs_ax + n_gs_cbar_gap + n_gs_cbar + n_gs_ax_gap)
        x1 = x0 + n_gs_ax
        ax = plt.subplot(gs[y0:y1, x0:x1])
        x0 = x1 + n_gs_cbar_gap
        x1 = x0 + n_gs_cbar
        cax = plt.subplot(gs[y0:y1, x0:x1])

        # Colour map
        cmap = plt.get_cmap(cmaps[i])
        norm = mpl.colors.Normalize(
            vmin=np.min(plot_vectors[i]), vmax=np.max(plot_vectors[i])
        )

        # Plot
        scatter = ax.scatter(
            A1_x,
            A1_y,
            c=plot_vectors[i],
            norm=norm,
            cmap=cmap,
            s=A1_size,
            edgecolors="none",
        )
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_facecolor((0.9, 0.9, 0.9))
        ax.set_xlim((0, boxsize_l))
        ax.set_ylim((0, boxsize_l))

        # Colour bar
        cbar = plt.colorbar(scatter, cax)
        cbar.ax.tick_params(labelsize=14)
        cbar.set_label(plot_quantity[i], rotation=90, labelpad=8, fontsize=18)

    plt.savefig("square_%s_%04d.png" % (type, snap), dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    # Load snapshot data
    mpl.use("Agg")
    type = sys.argv[1]
    snap = int(sys.argv[2])
    assert type in ["equal_spacing", "equal_mass"]
    snap_file = "square_%s_%04d.hdf5" % (type, snap)

    with h5py.File(snap_file, "r") as f:
        boxsize_l = f["Header"].attrs["BoxSize"][0]
        A1_x = f["/PartType0/Coordinates"][:, 0]
        A1_y = f["/PartType0/Coordinates"][:, 1]
        A1_z = f["/PartType0/Coordinates"][:, 2]
        A1_rho = f["/PartType0/Densities"][:]
        A1_u = f["/PartType0/InternalEnergies"][:]
        A1_P = f["/PartType0/Pressures"][:]
        A1_m = f["/PartType0/Masses"][:]

    # Sort arrays based on z position
    sort_indices = np.argsort(A1_z)
    A1_x = A1_x[sort_indices]
    A1_y = A1_y[sort_indices]
    A1_z = A1_z[sort_indices]
    A1_rho = A1_rho[sort_indices]
    A1_u = A1_u[sort_indices]
    A1_P = A1_P[sort_indices]
    A1_m = A1_m[sort_indices]

    # Mask to select slice
    slice_thickness = 0.1
    slice_pos_z = 0.5 * (np.max(A1_z) + np.min(A1_z))
    mask_slice = np.logical_and(
        A1_z > slice_pos_z - 0.5 * slice_thickness,
        A1_z < slice_pos_z + 0.5 * slice_thickness,
    )

    # Select particles to plot
    A1_x_slice = A1_x[mask_slice]
    A1_y_slice = A1_y[mask_slice]
    A1_rho_slice = A1_rho[mask_slice]
    A1_u_slice = A1_u[mask_slice]
    A1_P_slice = A1_P[mask_slice]
    A1_m_slice = A1_m[mask_slice]

    # Size of plotted particles
    size_factor = 5e4
    A1_size = size_factor * (A1_m_slice / A1_rho_slice) ** (2 / 3)

    # Plot figure
    plot_square(
        A1_x_slice, A1_y_slice, A1_rho_slice, A1_u_slice, A1_P_slice, A1_size, boxsize_l
    )
