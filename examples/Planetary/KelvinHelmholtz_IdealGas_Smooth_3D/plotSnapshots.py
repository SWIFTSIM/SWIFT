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
Generate plot of the 3D Kelvin--Helmholtz instability.
"""

import sys
import h5py
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


def make_axes():
    # Use Gridspec to set up figure
    n_gs_ax = 40
    n_gs_ax_gap = 1
    n_gs_cbar_gap = 1
    n_gs_cbar = 2

    ax_len = 5
    ax_gap_len = n_gs_ax_gap * ax_len / n_gs_ax
    cbar_gap_len = n_gs_cbar_gap * ax_len / n_gs_ax
    cbar_len = n_gs_cbar * ax_len / n_gs_ax

    fig = plt.figure(
        figsize=(3 * ax_len + 2 * ax_gap_len + cbar_gap_len + cbar_len, ax_len)
    )
    gs = mpl.gridspec.GridSpec(
        nrows=n_gs_ax, ncols=3 * n_gs_ax + 2 * n_gs_ax_gap + n_gs_cbar_gap + n_gs_cbar
    )

    ax_0 = plt.subplot(gs[:n_gs_ax, :n_gs_ax])
    ax_1 = plt.subplot(gs[:n_gs_ax, n_gs_ax + n_gs_ax_gap : 2 * n_gs_ax + n_gs_ax_gap])
    ax_2 = plt.subplot(
        gs[:n_gs_ax, 2 * n_gs_ax + 2 * n_gs_ax_gap : 3 * n_gs_ax + 2 * n_gs_ax_gap]
    )
    cax = plt.subplot(
        gs[
            :n_gs_ax,
            3 * n_gs_ax
            + 2 * n_gs_ax_gap
            + n_gs_cbar_gap : 3 * n_gs_ax
            + 2 * n_gs_ax_gap
            + n_gs_cbar_gap
            + n_gs_cbar,
        ]
    )

    axs = [ax_0, ax_1, ax_2]

    return axs, cax


def plot_kh(ax, snap, cmap, norm):

    # Load data
    snap_file = "kelvin_helmholtz_%04d.hdf5" % snap

    with h5py.File(snap_file, "r") as f:
        boxsize_l = f["Header"].attrs["BoxSize"][0]
        A1_x = f["/PartType0/Coordinates"][:, 0]
        A1_y = f["/PartType0/Coordinates"][:, 1]
        A1_z = f["/PartType0/Coordinates"][:, 2]
        A1_rho = f["/PartType0/Densities"][:]
        A1_m = f["/PartType0/Masses"][:]

    # Sort arrays based on z position
    sort_indices = np.argsort(A1_z)
    A1_x = A1_x[sort_indices]
    A1_y = A1_y[sort_indices]
    A1_z = A1_z[sort_indices]
    A1_rho = A1_rho[sort_indices]
    A1_m = A1_m[sort_indices]

    # Mask to select slice
    slice_thickness = 0.2 * (np.max(A1_z) - np.min(A1_z))
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

    # Size of plotted particles
    size_factor = 5e4
    A1_size = size_factor * (A1_m_slice / A1_rho_slice) ** (2 / 3) / boxsize_l ** 2

    # Plot
    scatter = ax.scatter(
        A1_x_slice,
        A1_y_slice,
        c=A1_rho_slice,
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


if __name__ == "__main__":

    # Set colormap
    cmap = plt.get_cmap("Spectral_r")
    vmin, vmax = 0.95, 2.05
    norm = mpl.colors.Normalize(vmin=vmin, vmax=vmax)

    # Generate axes
    axs, cax = make_axes()

    # The three snapshots to be plotted
    snaps = [15, 20, 25]
    times = ["1.5", "2.0", "2.5"]

    # Plot
    for i, snap in enumerate(snaps):
        ax = axs[i]
        time = times[i]

        plot_kh(ax, snap, cmap, norm)
        ax.text(
            0.5,
            -0.1,
            r"$t =\;$" + time + r"$\, \tau^{}_{\rm KH}$",
            horizontalalignment="center",
            size=18,
            transform=ax.transAxes,
        )

    # Colour bar
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    cbar = plt.colorbar(sm, cax)
    cbar.ax.tick_params(labelsize=14)
    cbar.set_label("Density", rotation=90, labelpad=8, fontsize=18)

    plt.savefig("kelvin_helmholtz.png", dpi=300, bbox_inches="tight")
