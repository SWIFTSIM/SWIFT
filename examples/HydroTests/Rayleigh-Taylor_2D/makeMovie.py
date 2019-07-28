###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
#                    Josh Borrow (joshua.borrow@durham.ac.uk)
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
from swiftsimio.visualisation import scatter
from matplotlib.animation import FuncAnimation
from matplotlib.colors import LogNorm
import matplotlib.pyplot as plt

from numpy import max, min
import numpy as np

try:
    # Try and load this, otherwise we're stuck with serial
    from p_tqdm import p_map
except:
    print("Try installing the p_tqdm package to make movie frames in parallel")
    p_map = map
    pass


info_frames = 40
generate_png = False
text_args = dict(color="black")


class Metadata(object):
    """
    Copy the useful data in order to decrease the memory usage
    """
    def __init__(self, data):
        metadata = data.metadata
        self.t = metadata.t
        try:
            self.viscosity_info = metadata.viscosity_info
        except:
            self.viscosity_info = "No info"
        try:
            self.diffusion_info = metadata.diffusion_info
        except:
            self.diffusion_info = "No info"

        self.code_info = metadata.code_info
        self.compiler_info = metadata.compiler_info
        self.hydro_info = metadata.hydro_info


def project(data, m_res, property, ylim):
    x, y, _ = data.gas.coordinates.value.T

    mask = np.logical_and(y >= ylim[0], y <= ylim[1])

    x = x[mask]
    y = y[mask] - np.float64(ylim[0])

    h = data.gas.smoothing_length[mask]

    if property == "density":
        property = "masses"

    if property is not None:
        quant = getattr(data.gas, property).value[mask]
    else:
        quant = np.ones_like(x)

    image = scatter(x=x, y=y, m=quant, h=h, res=m_res)

    return image.T


def load_and_make_image(filename, res, property):
    image = np.zeros(res, dtype=np.float32)
    m_res = min(res)
    border = int(0.2 * m_res)

    # first part of the image
    ylim = np.array([0., 1.])

    data = load(filename)
    image[:m_res, :m_res] = project(data, m_res, property, ylim)
    if property != "density":
        image[:m_res, :m_res] /= project(data, m_res, None, ylim)

    # second part of the image
    ylim = np.array([0.5, 1.5])

    left = -m_res + border
    image[left:, :] = project(data, m_res, property, ylim)[left:, :]
    if property != "density":
        image[left:, :] /= project(data, m_res, None, ylim)[left:, :]

    metadata = Metadata(data)
    return image, metadata


def time_formatter(metadata):
    return f"$t = {metadata.t:2.2f}$"


def create_movie(filename, start, stop, resolution, property, output_filename):
    """
    Creates the movie with:

    snapshots named {filename}_{start}.hdf5 to {filename}_{stop}.hdf5
    at {resolution} x {resolution} pixel size and smoothing the given
    {property} and outputting to {output_filename}.
    """

    if property != "density":
        name = property
    else:
        name = "Fluid Density $\\rho$"

    def baked_in_load(n):
        f = filename + "_{:04d}.hdf5".format(n)
        return load_and_make_image(f, resolution, property)

    # Make frames in parallel (reading also parallel!)
    frames, metadata = zip(*p_map(baked_in_load, list(range(start, stop))))

    vmax = max(list(p_map(max, frames)))
    vmin = min(list(p_map(min, frames)))

    fig, ax = plt.subplots(figsize=(8, 1.5 * 8), dpi=resolution[0] // 8)
    ax.axis("off")
    fig.subplots_adjust(0, 0, 1, 1)

    norm = LogNorm(vmin=vmin, vmax=vmax, clip="black")

    image = ax.imshow(np.zeros_like(frames[0]), origin="lower",
                      norm=norm)

    description_text = ax.text(
        0.5, 0.5,
        get_simulation_information(metadata[0]),
        va="center", ha="center",
        **text_args, transform=ax.transAxes,
    )

    time_text = ax.text(
        (1 - 0.025 * 0.25), 0.975,
        time_formatter(metadata[0]),
        **text_args,
        va="top", ha="right",
        transform=ax.transAxes,
    )

    ax.text(
        0.025 * 0.25, 0.975, name, **text_args, va="top", ha="left",
        transform=ax.transAxes
    )

    def frame(n):
        if n >= 0:
            image.set_array(frames[n])
            description_text.set_text("")
            time_text.set_text(time_formatter(metadata[n]))

        if generate_png:
            name = filename + "_{:04d}".format(n+info_frames)
            fig.savefig(name + ".png")

        else:
            return (image,)

    if generate_png:
        for i in range(-info_frames, stop-start):
            frame(i)
    else:
        animation = FuncAnimation(fig, frame,
                                  range(-info_frames, stop-start),
                                  interval=40)
        animation.save(output_filename)


def get_simulation_information(metadata):
    """
    Generates a string from the SWIFT metadata
    """
    viscosity = metadata.viscosity_info
    diffusion = metadata.diffusion_info

    output = (
        "$\\bf{Rayleigh-Taylor}$ $\\bf{Instability}$\n\n"
        "$\\bf{SWIFT}$\n"
        + metadata.code_info
        + "\n\n"
        + "$\\bf{Compiler}$\n"
        + metadata.compiler_info
        + "\n\n"
        + "$\\bf{Hydrodynamics}$\n"
        + metadata.hydro_info
        + "\n\n"
        + "$\\bf{Viscosity}$\n"
        + viscosity
        + "\n\n"
        + "$\\bf{Diffusion}$\n"
        + diffusion
    )

    return output


if __name__ == "__main__":
    from argparse import ArgumentParser

    parser = ArgumentParser(description="Creates a movie of the whole box.")

    parser.add_argument(
        "-s",
        "--snapshot",
        help="Snapshot name. Default: rayleigh_taylor",
        type=str,
        default="rayleigh_taylor",
    )

    parser.add_argument(
        "-i", "--initial", help="Initial snapshot. Default: 0",
        type=int, default=0
    )

    parser.add_argument(
        "-f", "--final", help="Final snapshot. Default: 40",
        type=int, default=40
    )

    parser.add_argument(
        "-o",
        "--output",
        help="Output filename. Default: rayleigh_taylor.mp4",
        type=str,
        default="rayleigh_taylor.mp4",
    )

    parser.add_argument(
        "-p",
        "--property",
        help="(swiftsimio) Property to plot. Default: density",
        type=str,
        default="density",
    )

    parser.add_argument(
        "-r",
        "--resolution",
        help="Resolution to make the movie at. Default: 512",
        type=int,
        default=512,
    )

    vars = parser.parse_args()

    yres = int(1.5 * vars.resolution)
    vars.resolution = [yres, vars.resolution]

    create_movie(
        filename=vars.snapshot,
        start=vars.initial,
        stop=vars.final,
        resolution=vars.resolution,
        property=vars.property,
        output_filename=vars.output,
    )
