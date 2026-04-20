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
Generate plot of the 3D fluid--solid interface.
"""

import h5py
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt


def load_snapshot(snap):
    # Load snapshot
    with h5py.File(f"fluid_solid_{snap:04d}.hdf5", "r") as f:
        coords = f["/PartType0/Coordinates"][:]
        vx = f["/PartType0/Velocities"][:, 0]
        rho = f["/PartType0/Densities"][:]
        m = f["/PartType0/Masses"][:]
        mat = f["/PartType0/MaterialIDs"][:]
    return coords, vx, rho, m, mat


def slice_data(coords, vx, rho, m, mat):
    # Extract a central slice in z for 2D visualisation
    z = coords[:, 2]
    z_mid = 0.5 * (z.min() + z.max())
    dz = 0.2 * (z.max() - z.min())
    mask = np.abs(z - z_mid) < 0.5 * dz
    return (
        coords[mask],
        vx[mask],
        rho[mask],
        m[mask],
        mat[mask],
    )


def plot(ax, coords, vx, rho, m, boxsize, cmap, norm):
    # Scatter plot of a 2D slice with particle-size scaling

    x = coords[:, 0]
    y = coords[:, 1]
    size = 5e4 * (m / rho) ** (2 / 3) / boxsize**2

    scatter = ax.scatter(
        x,
        y,
        c=vx,
        cmap=cmap,
        norm=norm,
        s=size,
        edgecolors="none",
    )
    ax.set_xlim(0, boxsize)
    ax.set_ylim(0, boxsize)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_facecolor((0.9, 0.9, 0.9))

    return scatter


def main():

    # Snapshots to plot (time evolution sequence)
    snaps = [0, 10, 15]
    times = ["0", "10e-5", "15e-5"]

    # Colormap and normalization for velocity field
    cmap = plt.get_cmap("Spectral_r")
    norm = mpl.colors.Normalize(vmin=0, vmax=1000)

    # Create panel layout
    fig, axs = plt.subplots(1, 3, figsize=(15, 5), constrained_layout=True)

    scatter_last = None

    for ax, snap, t in zip(axs, snaps, times):

        # Load and slice snapshot data
        coords, vx, rho, m, mat = load_snapshot(snap)
        coords, vx, rho, m, mat = slice_data(coords, vx, rho, m, mat)

        # Plot slice
        boxsize = coords[:, 0].max()
        scatter = plot(ax, coords, vx, rho, m, boxsize, cmap, norm)

        # Annotate time
        ax.text(
            0.5,
            -0.1,
            rf"$t = {t}\,$s",
            ha="center",
            transform=ax.transAxes,
            fontsize=14,
        )

        # Store last scatter for cbar
        scatter_last = scatter

    # Shared cbar across all panels
    cbar = fig.colorbar(scatter_last, ax=axs, shrink=0.8)
    cbar.set_label(r"$v_x$")

    # Save figure
    plt.savefig("fluid_solid.png", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()