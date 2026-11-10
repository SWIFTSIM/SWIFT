###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 SWIFT developers.
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
Makes a movie of the dam break: the density averaged through the thin (y)
direction of the tank, for every snapshot found in the current directory.

Requires swiftsimio and ffmpeg.
"""

from glob import glob

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import FuncAnimation
from swiftsimio import cosmo_array, load
from swiftsimio.visualisation.projection import project_gas

resolution = 1024
snapshot_name = "dam_break"
cmap = "Blues"
fps = 25

# View along the y axis: (x, y, z) -> (x, z, -y)
rotation_matrix = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]])


def get_image(filename):
    """
    Gets the y-averaged density image for a snapshot, and also returns the
    associated simulation time and number of particles.
    """
    data = load(filename)
    coords = data.gas.coordinates
    boxsize = data.metadata.boxsize

    def length(values):
        return cosmo_array(
            values,
            coords.units,
            comoving=coords.comoving,
            scale_factor=data.metadata.scale_factor,
            scale_exponent=1,
        )

    # Surface density of the tank seen along y, i.e. spanning x and z. Project
    # onto a square region and crop it, as only square images keep their units.
    size_x = boxsize[0].to_value(coords.units)
    size_z = boxsize[2].to_value(coords.units)
    size = max(size_x, size_z)
    image = project_gas(
        data,
        resolution,
        "masses",
        region=length([0.0, size, 0.0, size]),
        rotation_matrix=rotation_matrix,
        rotation_center=length([0.0, 0.0, 0.0]),
        periodic=False,
    )
    image = image[
        : int(np.ceil(resolution * size_x / size)),
        : int(np.ceil(resolution * size_z / size)),
    ]

    # Divide by the depth of the tank for the mean density through the slab.
    # The image is indexed [x, z]; imshow wants [z, x].
    density = image.to_value("kg/m**2").T / boxsize[1].to_value("m")
    extent = [0.0, boxsize[0].to_value("m"), 0.0, boxsize[2].to_value("m")]

    return density, data.metadata.t.to_value("s"), coords.shape[0], extent


if __name__ == "__main__":
    filenames = sorted(glob(f"{snapshot_name}_[0-9][0-9][0-9][0-9].hdf5"))
    if len(filenames) == 0:
        raise FileNotFoundError(f"No {snapshot_name}_XXXX.hdf5 snapshots found")

    images, times, counts, extents = zip(*[get_image(f) for f in filenames])

    fig, ax = plt.subplots(figsize=(8, 8 * images[0].shape[0] / images[0].shape[1]))
    fig.subplots_adjust(0.1, 0.1, 0.9, 0.95)
    image = ax.imshow(
        images[0], extent=extents[0], origin="lower", cmap=cmap, vmin=0.0, vmax=1100.0
    )
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z [m]")
    fig.colorbar(image, ax=ax, label=r"Mean density [kg m$^{-3}$]", shrink=0.8)
    title = ax.set_title("")

    def animate(n):
        image.set_array(images[n])
        title.set_text(f"t = {times[n]:.4f} s, {counts[n]} particles")
        return (image,)

    animation = FuncAnimation(fig, animate, range(len(images)), interval=1000 / fps)
    animation.save(filename=f"{snapshot_name}.mp4", fps=fps, dpi=150)
    print(f"Wrote {snapshot_name}.mp4 ({len(images)} frames)")
