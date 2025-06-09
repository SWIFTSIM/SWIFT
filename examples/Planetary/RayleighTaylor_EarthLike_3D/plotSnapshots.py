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
Generate plot of the 3D Rayleigh--Taylor instability.
"""

import sys
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


def make_axes():
    # Use Gridspec to set up figure
    n_gs_ax_x = 40
    n_gs_ax_y = 80
    n_gs_ax_gap = 1
    n_gs_cbar_gap = 1
    n_gs_cbar = 2

    ax_len_x = 5
    ax_len_y = 10
    ax_gap_len = n_gs_ax_gap * ax_len_x / n_gs_ax_x
    cbar_gap_len = n_gs_cbar_gap * ax_len_x / n_gs_ax_x
    cbar_len = n_gs_cbar * ax_len_x / n_gs_ax_x

    fig = plt.figure(
        figsize=(3 * ax_len_x + 2 * ax_gap_len + cbar_gap_len + cbar_len, ax_len_y)
    )
    gs = mpl.gridspec.GridSpec(
        nrows=n_gs_ax_y,
        ncols=3 * n_gs_ax_x + 2 * n_gs_ax_gap + n_gs_cbar_gap + n_gs_cbar,
    )

    ax_0 = plt.subplot(gs[:n_gs_ax_y, :n_gs_ax_x])
    ax_1 = plt.subplot(
        gs[:n_gs_ax_y, n_gs_ax_x + n_gs_ax_gap : 2 * n_gs_ax_x + n_gs_ax_gap]
    )
    ax_2 = plt.subplot(
        gs[
            :n_gs_ax_y,
            2 * n_gs_ax_x + 2 * n_gs_ax_gap : 3 * n_gs_ax_x + 2 * n_gs_ax_gap,
        ]
    )

    cax_0 = plt.subplot(
        gs[
            : int(0.5 * n_gs_ax_y) - 1,
            3 * n_gs_ax_x
            + 2 * n_gs_ax_gap
            + n_gs_cbar_gap : 3 * n_gs_ax_x
            + 2 * n_gs_ax_gap
            + n_gs_cbar_gap
            + n_gs_cbar,
        ]
    )
    cax_1 = plt.subplot(
        gs[
            int(0.5 * n_gs_ax_y) + 1 :,
            3 * n_gs_ax_x
            + 2 * n_gs_ax_gap
            + n_gs_cbar_gap : 3 * n_gs_ax_x
            + 2 * n_gs_ax_gap
            + n_gs_cbar_gap
            + n_gs_cbar,
        ]
    )

    axs = [ax_0, ax_1, ax_2]
    caxs = [cax_0, cax_1]

    return axs, caxs


def plot_kh(ax, snap, mat_id1, mat_id2, cmap1, cmap2, norm1, norm2):

    # Load data
    snap_file = "rayleigh_taylor_%04d.hdf5" % snap

    with h5py.File(snap_file, "r") as f:
        # Units from file metadata to SI
        m = float(f["Units"].attrs["Unit mass in cgs (U_M)"][0]) * 1e-3
        l = float(f["Units"].attrs["Unit length in cgs (U_L)"][0]) * 1e-2

        boxsize_x = f["Header"].attrs["BoxSize"][0] * l
        boxsize_y = f["Header"].attrs["BoxSize"][1] * l
        A1_x = f["/PartType0/Coordinates"][:, 0] * l
        A1_y = f["/PartType0/Coordinates"][:, 1] * l
        A1_z = f["/PartType0/Coordinates"][:, 2] * l
        A1_rho = f["/PartType0/Densities"][:] * (m / l ** 3)
        A1_m = f["/PartType0/Masses"][:] * m
        A1_mat_id = f["/PartType0/MaterialIDs"][:]

    # Sort arrays based on z position
    sort_indices = np.argsort(A1_z)
    A1_x = A1_x[sort_indices]
    A1_y = A1_y[sort_indices]
    A1_z = A1_z[sort_indices]
    A1_rho = A1_rho[sort_indices]
    A1_m = A1_m[sort_indices]
    A1_mat_id = A1_mat_id[sort_indices]

    # Mask to select slice
    slice_thickness = 0.1 * (np.max(A1_z) - np.min(A1_z))
    slice_pos_z = 0.5 * (np.max(A1_z) + np.min(A1_z))
    mask_slice = np.logical_and(
        A1_z > slice_pos_z - 0.5 * slice_thickness,
        A1_z < slice_pos_z + 0.5 * slice_thickness,
    )

    # Select particles to plot
    A1_x_slice = A1_x[mask_slice]
    A1_y_slice = A1_y[mask_slice]
    A1_rho_slice = A1_rho[mask_slice]
    A1_m_slice = A1_m[mask_slice]
    A1_mat_id_slice = A1_mat_id[mask_slice]

    # Size of plotted particles
    size_factor = 5e4
    A1_size = size_factor * (A1_m_slice / A1_rho_slice) ** (2 / 3) / boxsize_x ** 2

    mask_mat1 = A1_mat_id_slice == mat_id1
    mask_mat2 = A1_mat_id_slice == mat_id2

    # Plot
    scatter = ax.scatter(
        A1_x_slice[mask_mat1],
        A1_y_slice[mask_mat1],
        c=A1_rho_slice[mask_mat1],
        norm=norm1,
        cmap=cmap1,
        s=A1_size[mask_mat1],
        edgecolors="none",
    )

    scatter = ax.scatter(
        A1_x_slice[mask_mat2],
        A1_y_slice[mask_mat2],
        c=A1_rho_slice[mask_mat2],
        norm=norm2,
        cmap=cmap2,
        s=A1_size[mask_mat2],
        edgecolors="none",
    )

    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor((0.9, 0.9, 0.9))
    ax.set_xlim((0.0, boxsize_x))
    ax.set_ylim((0.05 * boxsize_y, 0.95 * boxsize_y))


if __name__ == "__main__":

    # Set colormap
    cmap1 = plt.get_cmap("YlOrRd")
    mat_id1 = 402
    rho_min1 = 7950
    rho_max1 = 10050
    norm1 = mpl.colors.Normalize(vmin=rho_min1, vmax=rho_max1)

    cmap2 = plt.get_cmap("Blues_r")
    mat_id2 = 400
    rho_min2 = 4950
    rho_max2 = 5550
    norm2 = mpl.colors.Normalize(vmin=rho_min2, vmax=rho_max2)

    # Generate axes
    axs, caxs = make_axes()

    # The three snapshots to be plotted
    snaps = [8, 12, 16]
    times = ["400", "600", "800"]

    # Plot
    for i, snap in enumerate(snaps):
        ax = axs[i]
        time = times[i]

        plot_kh(ax, snap, mat_id1, mat_id2, cmap1, cmap2, norm1, norm2)
        ax.text(
            0.5,
            -0.05,
            r"$t =\;$" + time + r"$\,$s",
            horizontalalignment="center",
            size=18,
            transform=ax.transAxes,
        )

    # Colour bar
    sm1 = plt.cm.ScalarMappable(cmap=cmap1, norm=norm1)
    cbar1 = plt.colorbar(sm1, caxs[0])
    cbar1.ax.tick_params(labelsize=14)
    cbar1.set_label(r"Iron density (kg/m$^3$)", rotation=90, labelpad=8, fontsize=12)

    sm2 = plt.cm.ScalarMappable(cmap=cmap2, norm=norm2)
    cbar2 = plt.colorbar(sm2, caxs[1])
    cbar2.ax.tick_params(labelsize=14)
    cbar2.set_label(r"Rock density (kg/m$^3$)", rotation=90, labelpad=16, fontsize=12)

    plt.savefig("rayleigh_taylor.png", dpi=300, bbox_inches="tight")
