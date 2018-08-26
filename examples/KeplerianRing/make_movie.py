"""
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2017
#
# Josh Borrow (joshua.borrow@durham.ac.uk)
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
#
# -----------------------------------------------------------------------------
#
# This program creates test plots for the initial condition generator provided
# for the Keplerian Ring example.
#
###############################################################################
"""

import matplotlib
matplotlib.use("Agg")


import matplotlib.pyplot as plt
import matplotlib.animation as anim
import numpy as np
import h5py as h5

from tqdm import tqdm


def get_colours(numbers, norm=None, cm_name="viridis"):
    if norm is None:
        norm = matplotlib.colors.Normalize(vmax=numbers.max(), vmin=numbers.min())

    cmap = matplotlib.cm.get_cmap(cm_name)

    return cmap(norm(numbers)), norm
        

def load_data(filename, silent=True, extra=None, logextra=True, exclude=None):
    if not silent:
        print(f"Loading data from {filename}")

    with h5.File(filename, "r") as file_handle:
        coords = file_handle['PartType0']['Coordinates'][...]
        time = float(file_handle['Header'].attrs['Time'])
        boxsize = file_handle['Header'].attrs['BoxSize']

        # For some old runs we have z=0
        if np.sum(coords[:, 2]) == 0:
            centre_of_box = np.array(list(boxsize[:2]/2) + [0])
        else:
            centre_of_box = boxsize/2

        
        if exclude is not None:
            distance_from_centre = coords - centre_of_box
            r2 = np.sum(distance_from_centre * distance_from_centre, 1)
            mask = r2 < (exclude * exclude)

            masked = np.ma.array(coords, mask=np.array([mask]*3))

            coords = masked.compressed()
        

        if extra is not None:
            other = file_handle['PartType0'][extra][...]
            if exclude is not None:
                masked = np.ma.array(other, mask=mask)

                other = masked.compressed()

            if logextra:
                other = np.log(other)

            return coords, time, other
        else:
            return coords, time


def rms(x):
    return np.sqrt(sum(x**2))


def rotation_velocity_at_r(r, params):
    """
    Gets the rotation velocity at a given radius r by assuming it is keplerian.

    Assumes we are in cgs units, which may one day be our downfall.
    """

    unit_length = float(params[r"InternalUnitSystem:UnitLength_in_cgs"])

    if unit_length != 1.:
        print(f"Your unit length: {unit_length}")
        raise InternalUnitSystemError(
            "This function is only created to handle CGS units."
        )

    central_mass = float(params["PointMassPotential:mass"])
    G = 6.67408e-8

    v = np.sqrt( G * central_mass / r)

    return v


def get_rotation_period_at_r(r, params):
    """
    Gets the rotation period at a given radius r, assuming a keplerian
    orbit.
    """
    v = rotation_velocity_at_r(r, params)

    return 2*np.pi / v


def get_metadata(filename, r=1):
    """ The metadata should be extracted from the first snapshot. """
    with h5.File(filename, "r") as file_handle:
        header = file_handle['Header'].attrs
        code = file_handle['Code'].attrs
        hydro = file_handle['HydroScheme'].attrs
        params = file_handle['Parameters'].attrs

        period = get_rotation_period_at_r(r, params)

        return_values = {
            "header" : dict(header),
            "code" : dict(code),
            "period" : float(period),
            "hydro" : dict(hydro),
            "params" : dict(params)
        }

    return return_values


def plot_single(number, scatter, text, metadata, ax, extra=None, norm=None):
    filename = "keplerian_ring_{:04d}.hdf5".format(number)

    if extra is not None:
        coordinates, time, other = load_data(filename, extra=extra)
    else:
        coordinates, time = load_data(filename)


    text.set_text(
        "Time: {:1.2f} | Rotations {:1.2f}".format(
            time,
            time/metadata['period'],
        )
    )

    data = coordinates[:, 0:2]
    scatter.set_offsets(data)

    if extra is not None:
        colours, _ = get_colours(other, norm)
        scatter.set_color(colours)

    return scatter,


if __name__ == "__main__":
    import os
    import argparse as ap

    parser = ap.ArgumentParser(
        description="""
                   Plots a movie of the current snapshots in your
                   directory. Can also colourmap your information.
                   """
    )

    parser.add_argument(
        "-e",
        "--exclude_central",
        help="""
             Region from the centre of the ring to exclude from the
             vmin/vmax of the colourmap. Note that these particles are
             still plotted -- they are just excluded for the purposes
             of colourmapping. Default: 1 simulation unit.
             """,
        default=1.,
        required=False
    )

    parser.add_argument(
        "-c",
        "--cmap",
        help="""
             The item from the GADGET hdf5 file to clourmap with.
             Examples include Density, InternalEnergy.
             Default: don't use a colourmap. (Much faster).
             """,
        required=False,
        default=None
    )

    parser.add_argument(
        "-m",
        "--max",
        help="""
             Maximum radii to plot.
             Default: 2.5.
             """,
        required=False,
        default=2.5,
    )
    
    args = vars(parser.parse_args())

    # Look for the number of files in the directory.
    i = 0
    while True:
        if os.path.isfile("keplerian_ring_{:04d}.hdf5".format(i)):
            i += 1
        else:
            break

        if i > 10000:
            break


   
    # Now deal with the colourmapping (if present)
    if args["cmap"] is not None:
        _, _, numbers0 = load_data("keplerian_ring_0000.hdf5", extra=args["cmap"], exclude=float(args["exclude_central"]))
        _, _, numbersend = load_data("keplerian_ring_{:04d}.hdf5".format(i-1), extra=args["cmap"], exclude=float(args["exclude_central"]))
        vmax = max([numbers0.max(), numbersend.max()])
        vmin = min([numbers0.min(), numbersend.min()])

        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    else:
        norm = None


    # Initial plot setup

    metadata = get_metadata("keplerian_ring_0000.hdf5")
    n_particle = metadata['header']['NumPart_Total'][0]

    fig = plt.figure(figsize=(8, 8))
    ax = fig.add_subplot(111)

    scatter = ax.scatter([0]*n_particle, [0]*n_particle, s=0.5, marker="o")
    diff = float(args["max"])
    left = metadata['header']['BoxSize'][0]/2 - diff
    right = metadata['header']['BoxSize'][0]/2 + diff
    ax.set_xlim(left, right)
    ax.set_ylim(left, right)

    offset = 0.25
    time_text = ax.text(
        offset + left,
        offset + left,
        "Time: {:1.2f} | Rotations {:1.2f}".format(
            0,
            0/metadata['period'],
        )
    )

    ax.text(
        offset + left,
        right-offset-0.35,
        "Code: {} {} | {} {} \nHydro {}\n$\eta$={:1.4f}".format(
            metadata['code']['Git Branch'].decode("utf-8"),
            metadata['code']['Git Revision'].decode("utf-8"),
            metadata['code']['Compiler Name'].decode("utf-8"),
            metadata['code']['Compiler Version'].decode("utf-8"),
            metadata['hydro']['Scheme'].decode("utf-8"),
            metadata['hydro']['Kernel eta'][0],
        )
    )

    ax.set_title("Keplerian Ring Test")
    ax.set_xlabel("$x$ position")
    ax.set_ylabel("$y$ position")

    
    anim = anim.FuncAnimation(
        fig,
        plot_single,
        tqdm(np.arange(i)),
        fargs = [
            scatter,
            time_text,
            metadata,
            ax,
            args["cmap"],
            norm
        ],
        interval=50,
        repeat_delay=3000,
        blit=True,
    )

    anim.save("keplerian_ring.mp4", dpi=int(640/8))
