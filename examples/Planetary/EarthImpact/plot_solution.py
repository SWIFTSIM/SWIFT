###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2019 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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

# Plot the snapshots from the example giant impact on the proto-Earth, showing
# the particles in a thin slice near z=0, coloured by their material.

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import swiftsimio as sw
import unyt

font_size = 20
params = {
    "axes.labelsize": font_size,
    "font.size": font_size,
    "xtick.labelsize": font_size,
    "ytick.labelsize": font_size,
    "font.family": "serif",
}
matplotlib.rcParams.update(params)

# Material IDs ( = type_id * type_factor + unit_id )
type_factor = 100
type_Til = 1
id_body = 10000
# Name and ID
Di_mat_id = {
    "Til_iron": type_Til * type_factor,
    "Til_iron_2": type_Til * type_factor + id_body,
    "Til_granite": type_Til * type_factor + 1,
    "Til_granite_2": type_Til * type_factor + 1 + id_body,
}
# Colour
Di_mat_colour = {
    "Til_iron": "darkgray",
    "Til_granite": "orangered",
    "Til_iron_2": "saddlebrown",
    "Til_granite_2": "gold",
}
Di_id_colour = {Di_mat_id[mat]: colour for mat, colour in Di_mat_colour.items()}


def load_snapshot(snapshot_id, ax_lim):
    """ Select and load the particles to plot. """
    # Snapshot to load
    snapshot = "earth_impact_%04d.hdf5" % snapshot_id

    # Only load data with the axis limits and below z=0
    ax_lim = 0.1
    mask = sw.mask(snapshot)
    box_mid = 0.5 * mask.metadata.boxsize[0].to(unyt.Rearth)
    x_min = box_mid - ax_lim * unyt.Rearth
    x_max = box_mid + ax_lim * unyt.Rearth
    load_region = [[x_min, x_max], [x_min, x_max], [x_min, box_mid]]
    mask.constrain_spatial(load_region)

    # Load
    data = sw.load(snapshot, mask=mask)
    pos = data.gas.coordinates.to(unyt.Rearth) - box_mid
    id = data.gas.particle_ids
    mat_id = data.gas.material_ids.value

    # Restrict to z < 0
    sel = np.where(pos[:, 2] < 0)[0]
    pos = pos[sel]
    id = id[sel]
    mat_id = mat_id[sel]

    # Sort in z order so higher particles are plotted on top
    sort = np.argsort(pos[:, 2])
    pos = pos[sort]
    id = id[sort]
    mat_id = mat_id[sort]

    # Edit material IDs for particles in the impactor
    num_in_target = 99740
    mat_id[num_in_target <= id] += id_body

    return pos, mat_id


def plot_snapshot(pos, mat_id, ax_lim):
    """ Plot the particles, coloured by their material. """
    plt.figure(figsize=(7, 7))
    ax = plt.gca()
    ax.set_aspect("equal")

    colour = np.empty(len(pos), dtype=object)
    for id_c, c in Di_id_colour.items():
        colour[mat_id == id_c] = c

    ax.scatter(
        pos[:, 0],
        pos[:, 1],
        c=colour,
        edgecolors="none",
        marker=".",
        s=10,
        alpha=0.5,
    )

    ax.set_xlim(-ax_lim, ax_lim)
    ax.set_yticks(ax.get_xticks())
    ax.set_ylim(-ax_lim, ax_lim)
    ax.set_xlabel(r"x Position ($R_\oplus$)")
    ax.set_ylabel(r"y Position ($R_\oplus$)")

    plt.tight_layout()


if __name__ == "__main__":
    print()
    # Axis limits (Earth radii)
    ax_lim = 3.4

    # Plot each snapshot
    for snapshot_id in range(37):
        # Load the data
        pos, mat_id = load_snapshot(snapshot_id, ax_lim)

        # Plot the data
        plot_snapshot(pos, mat_id, ax_lim)

        # Save the figure
        save = "earth_impact_%04d.png" % snapshot_id
        plt.savefig(save, dpi=100)

        print("\rSaved %s" % save)

    print()
