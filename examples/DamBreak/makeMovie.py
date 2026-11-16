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
Makes a movie of the dam break seen along the thin (y) direction of the tank,
for every snapshot found in the current directory.

By default the movie shows the density averaged through the tank, which
includes any empty space. With --quantity density or --quantity pressure it
instead shows the mass-weighted mean of that quantity over the particles along
each line of sight. The colour scale of the density by default runs from the
smallest to the largest density of any particle in any snapshot. That of the
pressure stops at the 99th percentile of the particle pressures, as a few
particles with extreme pressures would otherwise wash out the rest.

The frames are drawn in parallel, one snapshot per process, and then encoded
into the movie by ffmpeg.

Requires swiftsimio and ffmpeg.
"""

import argparse
import os
import subprocess
import tempfile
from functools import partial
from glob import glob
from multiprocessing import Pool

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from swiftsimio import cosmo_array, load
from swiftsimio.visualisation.projection import project_gas

resolution = 1024
snapshot_name = "dam_break"
fps = 25

# For each quantity that can be shown: the particle field it is made from
# (None for the density averaged through the tank), its units, the label of
# the colour bar, the colour map, the name of the movie, and the percentile
# of the particle values at which the colour scale ends by default.
quantities = {
    "mean_density": (
        None,
        "kg/m**3",
        r"Mean density [kg m$^{-3}$]",
        "Blues",
        f"{snapshot_name}.mp4",
        None,
    ),
    "density": (
        "densities",
        "kg/m**3",
        r"Density [kg m$^{-3}$]",
        "viridis",
        f"{snapshot_name}_density.mp4",
        100.0,
    ),
    "pressure": (
        "pressures",
        "Pa",
        "Pressure [Pa]",
        "viridis",
        f"{snapshot_name}_pressure.mp4",
        99.0,
    ),
}

# View along the y axis: (x, y, z) -> (x, z, -y)
rotation_matrix = np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 1.0], [0.0, -1.0, 0.0]])


def get_limits(filename, quantity):
    """
    Gets the default limits of the colour scale of a quantity for a snapshot:
    its smallest value over the particles, and its largest or a percentile.
    """
    field, units = quantities[quantity][:2]
    values = getattr(load(filename).gas, field).to_value(units)
    return values.min(), np.percentile(values, quantities[quantity][5])


def get_image(filename, quantity):
    """
    Gets the image of a quantity for a snapshot, and also returns the
    associated simulation time and number of particles.
    """
    field, units = quantities[quantity][:2]
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

    # Project along y, i.e. spanning x and z, onto a square region and crop
    # it, as only square images keep their units.
    size_x = boxsize[0].to_value(coords.units)
    size_z = boxsize[2].to_value(coords.units)
    size = max(size_x, size_z)

    def project(name):
        image = project_gas(
            data,
            resolution,
            name,
            region=length([0.0, size, 0.0, size]),
            rotation_matrix=rotation_matrix,
            rotation_center=length([0.0, 0.0, 0.0]),
            periodic=False,
        )
        return image[
            : int(np.ceil(resolution * size_x / size)),
            : int(np.ceil(resolution * size_z / size)),
        ]

    mass = project("masses")
    if field is None:
        # Divide the surface density by the depth of the tank for the mean
        # density through the slab.
        image = mass.to_value("kg/m**2") / boxsize[1].to_value("m")
    else:
        # Mass-weighted mean along the line of sight, left blank where there
        # are no particles.
        data.gas.weighted = data.gas.masses * getattr(data.gas, field)
        with np.errstate(invalid="ignore", divide="ignore"):
            image = (project("weighted") / mass).to_value(units)

    # The image is indexed [x, z]; imshow wants [z, x].
    extent = [0.0, boxsize[0].to_value("m"), 0.0, boxsize[2].to_value("m")]

    return image.T, data.metadata.t.to_value("s"), coords.shape[0], extent


def draw_frame(job, quantity, limits):
    """
    Draws the frame of the movie for one snapshot and saves it as an image.
    """
    filename, frame_name = job
    label, cmap = quantities[quantity][2:4]
    values, time, count, extent = get_image(filename, quantity)

    fig, ax = plt.subplots(figsize=(8, 8 * values.shape[0] / values.shape[1]))
    fig.subplots_adjust(0.1, 0.1, 0.9, 0.95)
    image = ax.imshow(
        values,
        extent=extent,
        origin="lower",
        cmap=cmap,
        vmin=limits[0],
        vmax=limits[1],
    )
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z [m]")
    fig.colorbar(image, ax=ax, label=label, shrink=0.8)
    ax.set_title(f"t = {time:.4f} s, {count} particles")
    fig.savefig(frame_name, dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    if hasattr(os, "sched_getaffinity"):
        available = len(os.sched_getaffinity(0))
    else:
        available = os.cpu_count()

    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--quantity",
        choices=quantities.keys(),
        default="mean_density",
        help="What to show (default: mean_density)",
    )
    parser.add_argument(
        "--limits",
        type=float,
        nargs=2,
        metavar=("MIN", "MAX"),
        help="Limits of the colour scale (default: 0 to 1100 for mean_density, "
        "smallest to largest particle value for density, smallest to 99th "
        "percentile for pressure)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=available,
        help="Number of snapshots to draw at once (default: all available cores)",
    )
    args = parser.parse_args()
    movie_name = quantities[args.quantity][4]

    # Check that ffmpeg runs now rather than after drawing every frame. It can
    # be on the path and still fail to start, e.g. for want of a library.
    try:
        subprocess.run(["ffmpeg", "-version"], check=True, capture_output=True)
    except (OSError, subprocess.CalledProcessError) as error:
        raise RuntimeError(
            "ffmpeg, which is needed to encode the movie, does not run"
        ) from error

    filenames = sorted(glob(f"{snapshot_name}_[0-9][0-9][0-9][0-9].hdf5"))
    if len(filenames) == 0:
        raise FileNotFoundError(f"No {snapshot_name}_XXXX.hdf5 snapshots found")

    with Pool(args.threads) as pool, tempfile.TemporaryDirectory(
        dir=".", prefix="frames_"
    ) as frame_dir:
        # The same colour scale is used for every frame.
        if args.limits is not None:
            limits = args.limits
        elif args.quantity == "mean_density":
            limits = (0.0, 1100.0)
        else:
            lows, highs = zip(
                *pool.map(partial(get_limits, quantity=args.quantity), filenames)
            )
            limits = (min(lows), max(highs))
            print(f"Colour scale from {limits[0]:.6g} to {limits[1]:.6g}")

        jobs = [
            (filename, os.path.join(frame_dir, f"frame_{n:04d}.png"))
            for n, filename in enumerate(filenames)
        ]
        draw = partial(draw_frame, quantity=args.quantity, limits=limits)
        for n, _ in enumerate(pool.imap_unordered(draw, jobs), start=1):
            print(f"\rDrew frame {n} of {len(jobs)}", end="", flush=True)
        print()

        # The pad filter rounds the frame size up to the even numbers of
        # pixels that the widely supported yuv420p pixel format requires.
        subprocess.run(
            [
                "ffmpeg",
                "-loglevel",
                "error",
                "-y",
                "-framerate",
                str(fps),
                "-i",
                os.path.join(frame_dir, "frame_%04d.png"),
                "-vf",
                "pad=ceil(iw/2)*2:ceil(ih/2)*2",
                "-pix_fmt",
                "yuv420p",
                movie_name,
            ],
            check=True,
        )
    print(f"Wrote {movie_name} ({len(filenames)} frames)")
