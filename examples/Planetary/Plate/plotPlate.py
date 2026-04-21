###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

import h5py
import numpy as np
import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable


for snap in range(400):
    # Load snapshot
    filename = "plate_%04d.hdf5" % snap

    with h5py.File(filename, "r") as sim:
        coords = sim["/PartType0/Coordinates"][:]
        vx = sim["/PartType0/Velocities"][:, 0]
        vy = sim["/PartType0/Velocities"][:, 1]
        P = sim["/PartType0/Pressures"][:]
        rho = sim["/PartType0/Densities"][:]
        ids = sim["/PartType0/ParticleIDs"][:]

    # Choose what to plot
    func_plotted = P

    # Plot
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    size = 5
    norm = mpl.colors.Normalize(vmin=-3e4, vmax=3e4)
    sm = plt.cm.ScalarMappable(cmap="rainbow", norm=norm)
    sm.set_array([])
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(sm, cax=cax, orientation="vertical")
    cbar.set_label("Pressure")

    ax.scatter(
        coords[:, 0],
        coords[:, 1],
        c=func_plotted,
        norm=norm,
        cmap="rainbow",
        alpha=1,
        s=size,
    )

    ax.set_aspect("equal", "box")
    ax.set_facecolor((0.4, 0.4, 0.4))
    ax.set_xlim((0, 0.3))
    ax.set_ylim((0.1, 0.3))
    ax.set_xticklabels([])
    ax.set_yticklabels([])

    # Save
    plt.savefig(
        "images/plate_%04d.png" % snap,
        dpi=300,
        bbox_inches="tight",
    )

    plt.close(fig)
