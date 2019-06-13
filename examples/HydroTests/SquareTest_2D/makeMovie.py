###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk)
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

from swiftsimio import load
from swiftsimio.visualisation import project_gas_pixel_grid
from matplotlib.animation import FuncAnimation
import matplotlib.pyplot as plt

from numpy import max, min

try:
    # Try and load this, otherwise we're stuck with serial
    from p_tqdm import p_map

    map = p_map
except:
    print("Try installing the p_tqdm package to make movie frames in parallel")
    pass


def load_and_make_image(number, res, property):
    filename = "square_{:04d}.hdf5".format(number)
    data = load(filename)
    image = project_gas_pixel_grid(data, res, property)
    image_none = project_gas_pixel_grid(data, res, None)

    return image / image_none


def create_movie(filename, start, stop, resolution, property, output_filename):
    """
    Creates the movie with:

    snapshots named {filename}_{start}.hdf5 to {filename}_{stop}.hdf5
    at {resolution} x {resolution} pixel size and smoothing the given
    {property} and outputting to {output_filename}.
    """

    def baked_in_load(n):
        return load_and_make_image(n, resolution, property)

    # Make frames in parallel (reading also parallel!)
    frames = map(baked_in_load, list(range(start, stop)))

    vmax = max(list(map(max, frames)))
    vmin = min(list(map(min, frames)))

    fig, ax = plt.subplots(figsize=(1, 1))
    ax.axis("off")
    fig.subplots_adjust(0, 0, 1, 1)

    image = ax.imshow(frames[0], origin="lower", vmax=vmax, vmin=vmin)

    def frame(n):
        image.set_array(frames[n])

        return (image,)

    animation = FuncAnimation(fig, frame, range(0, 40), interval=40)
    animation.save(output_filename, dpi=resolution)


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="""Creates a movie of the whole box.""")

    parser.add_argument(
        "-s",
        "--snapshot",
        help="Snapshot name. Default: square",
        type=str,
        default="square",
    )

    parser.add_argument(
        "-i", "--initial", help="Initial snapshot. Default: 0", type=int, default=0
    )

    parser.add_argument(
        "-f", "--final", help="Final snapshot. Default: 40", type=int, default=40
    )

    parser.add_argument(
        "-o",
        "--output",
        help="Output filename. Default: EnergyMovie.mp4",
        type=str,
        default="EnergyMovie.mp4",
    )

    parser.add_argument(
        "-p",
        "--property",
        help="(swiftsimio) Property to plot. Default: internal_energy",
        type=str,
        default="internal_energy",
    )

    parser.add_argument(
        "-r",
        "--resolution",
        help="Resolution to make the movie at. Default: 512",
        type=int,
        default=512,
    )

    vars = parser.parse_args()

    create_movie(
        filename=vars.snapshot,
        start=vars.initial,
        stop=vars.final,
        resolution=vars.resolution,
        property=vars.property,
        output_filename=vars.output,
    )

    exit(0)

