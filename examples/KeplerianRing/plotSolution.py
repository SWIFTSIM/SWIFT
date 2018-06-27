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
###############################################################################
"""

# Plotting script for the Keplerian Ring example.
# We use yt for the projection rendering of the ring,
# and then our own density as a function of radius calculation.

import matplotlib
matplotlib.use("Agg")
matplotlib.rc("text", usetex=True)

try:
    import yt
    ytavail = True
except ImportError:
    print("yt not found. Falling back on homebrew plots.")
    ytavail = False

import h5py
import os

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec

try:
    from tqdm import tqdm
except ImportError:
    tqdm = lambda x: x

from make_movie import get_metadata


yt.funcs.mylog.setLevel(50)
tqdm.monitor_interval = 0


class InternalUnitSystemError(Exception):
    pass


def get_axes_grid(figure):
    """
    Grab our axes grid.
    """
    gs = gridspec.GridSpec(2, 30)

    grid = []
    
    grid.append(figure.add_subplot(gs[0, 0:10]))
    grid.append(figure.add_subplot(gs[0, 10:20]))
    grid.append(figure.add_subplot(gs[0, 20:]))
    grid.append(figure.add_subplot(gs[1, 0:9]))
    grid.append(figure.add_subplot(gs[1, 11:20]))
    grid.append(figure.add_subplot(gs[1, 21:]))

    return grid


def get_yt_actual_data(plot, name="density"):
    """
    Extracts the image data and colourmap from a yt plot.
    
    This is used to put on our own grid.
    """

    data = plot.plots[name].image.get_array()
    cmap = plot.plots[name].image.cmap

    return data, cmap


def chi_square(observed, expected):
    """
    The chi squared statistic.
        
    This also looks for where expected == 0 and masks over those particles to
    avoid divide by zero errors and unrealistically high chi squared.
    """

    mask = np.array(expected) != 0

    masked_expected = np.array(expected)[mask]
    masked_observed = np.array(observed)[mask]

    return sum(((masked_observed - masked_expected)**2)/masked_expected**2)


def load_data(filename):
    """
    Loads the data and extracts the relevant information for
    calculating the chi squared statistic and density(r) profiles.
    """

    with h5py.File(filename, "r") as file:
        boxsize = np.array(file["Header"].attrs["BoxSize"])

        # Check if z = 0 for all particles. If so we need to set the cetnre
        # to have z = 0 also.
        if np.sum(file["PartType0"]["Coordinates"][:, 2]) == 0:
            centre = [boxsize[0] / 2., boxsize[0] / 2., 0.]
        else:
            centre = boxsize / 2.

        radii = np.sqrt(np.sum(((file["PartType0"]["Coordinates"][...] - centre).T)**2, 0))
        masses = file["PartType0"]["Masses"][...]

    return radii, masses


def bin_density_r(radii, density, binrange, binnumber):
    """
    Bins the density as a funciton of radius.
    """

    bins = np.linspace(*binrange, binnumber)
    indicies = np.digitize(radii, bins)

    binned_masses = np.zeros(len(bins) - 1)
    
    for index, bin in enumerate(indicies):
        if bin >= len(bins) - 1:
            continue

        binned_masses[bin] += density[index]

    areas = [np.pi * (a**2 - b**2) for a, b in zip(bins[1:], bins[:-1])]
    binned_densities = binned_masses/areas
    
    return bins, binned_densities


def get_density_r(snapshot, filename="keplerian_ring", binrange=(0, 5), binnumber=50):
    """
    Gets the binned density as a function of radius.
    """
    snap = "{:04d}".format(snapshot)
    filename = f"{filename}_{snap}.hdf5"

    data = load_data(filename)

    return bin_density_r(*data, binrange, binnumber)


def get_mass_outside_inside(snap, radius=2, filename="keplerian_ring"):
    """
    Finds the mass outside and inside a radius. This is used to look at mass
    flow inside and outside of the ring with a postprocessing routine in
    get_derived_data.
    """
    snapshot = "{:04d}".format(snap)
    filename = f"{filename}_{snapshot}.hdf5"

    radii, masses = load_data(filename)

    below = sum(masses[radii < radius])
    above = sum(masses[radii > radius])

    return below, above


def get_derived_data(minsnap, maxsnap, filename="keplerian_ring", binnumber=50, radius=2):
    """
    Gets the derived data from our snapshots, i.e. the
    density(r) profile and the chi squared (based on the
    difference between the minsnap and the current snapshot).
    """

    initial = get_density_r(minsnap, filename, binnumber=binnumber)
    other_densities = [
        get_density_r(snap, binnumber=binnumber)[1] for snap in tqdm(
            range(minsnap+1, maxsnap+1), desc="Densities"
        )
    ]
    densities = [initial[1]] + other_densities

    masses_inside_outside = [
        get_mass_outside_inside(snap, radius=radius, filename=filename)[0] for snap in tqdm(
            range(minsnap, maxsnap+1), desc="Mass Flow"
        )
    ]

    # Between the initial conditions and the first snapshot we hope that there
    # has been no mass flow, hence the [0.] + 
    mass_flows = [0.] + [
        y - x for x, y in zip(
            masses_inside_outside[1:],
            masses_inside_outside[:-1]
        )
    ]

    cumulative_mass_flows = [
        sum(mass_flows[:x])/masses_inside_outside[0] for x in range(len(mass_flows))
    ]

    chisq = [
        chi_square(
            dens[int(0.1*binnumber):int(0.4*binnumber)],
            initial[1][int(0.1*binnumber):int(0.4*binnumber)]
        ) for dens in tqdm(densities, desc="Chi Squared")
    ]

    return initial[0], densities, chisq, cumulative_mass_flows


def plot_chisq(ax, minsnap, maxsnap, chisq, filename="keplerian_ring"):
    """
    Plot the chisq(rotation).
    """
    snapshots = np.arange(minsnap, maxsnap + 1)
    rotations = [convert_snapshot_number_to_rotations_at(1, snap, filename) for snap in snapshots]
    ax.plot(rotations, np.array(chisq)/max(chisq))

    ax.set_xlabel("Number of rotations")
    ax.set_ylabel("$\chi^2 / \chi^2_{{max}}$ = {:3.5f}".format(max(chisq)))

    return


def plot_mass_flow(ax, minsnap, maxsnap, mass_flow, filename="keplerian_ring"):
    """
    Plot the mass_flow(rotation).
    """
    snapshots = np.arange(minsnap, maxsnap + 1)
    rotations = [convert_snapshot_number_to_rotations_at(1, snap, filename) for snap in snapshots]

    ax.plot(rotations, mass_flow)

    ax.set_xlabel("Number of rotations")
    ax.set_ylabel(r"Mass flow out of ring ($M_{\rm ring}$)")

    return


def plot_density_r(ax, bins, densities, snaplist, filename="keplerian_ring"):
    """
    Make the density(r) plots.

    Densities is the _full_ list of density profiles, and
    snaplist is the ones that you wish to plot.
    """
    radii = [(x + y)/2 for x, y in zip(bins[1:], bins[:-1])]

    for snap in snaplist:
        index = snap - snaplist[0]
        rotations = convert_snapshot_number_to_rotations_at(1, snap, filename)
        ax.plot(radii, densities[index], label="{:2.2f} Rotations".format(rotations))

    ax.legend()
    ax.set_xlabel("Radius")
    ax.set_ylabel("Azimuthally Averaged Surface Density")

    return


def plot_extra_info(ax, filename):
    """
    Plots all of the extra information on the final axis.

    Give it a filename of any of the snapshots.
    """

    metadata = get_metadata(filename)
    
    git = metadata['code']['Git Revision'].decode("utf-8")
    compiler_name = metadata['code']['Compiler Name'].decode("utf-8")
    compiler_version = metadata['code']['Compiler Version'].decode("utf-8")
    scheme = metadata['hydro']['Scheme'].decode("utf-8")
    kernel = metadata['hydro']['Kernel function'].decode("utf-8")
    gas_gamma = metadata["hydro"]["Adiabatic index"][0]
    neighbors = metadata["hydro"]["Kernel target N_ngb"][0]
    eta = metadata["hydro"]["Kernel eta"][0]
    

    ax.text(-0.49, 0.9, "Keplerian Ring with  $\\gamma={:4.4f}$ in 2/3D".format(gas_gamma), fontsize=11)
    ax.text(-0.49, 0.8, f"Compiler: {compiler_name} {compiler_version}", fontsize=10)
    ax.text(-0.49, 0.7, "Rotations are quoted at $r=1$", fontsize=10)
    ax.plot([-0.49, 0.1], [0.62, 0.62], 'k-', lw=1)
    ax.text(-0.49, 0.5, f"$\\textsc{{Swift}}$ {git}", fontsize=10)
    ax.text(-0.49, 0.4, scheme, fontsize=10)
    ax.text(-0.49, 0.3, kernel, fontsize=10)
    ax.text(-0.49, 0.2, "${:2.2f}$ neighbours ($\\eta={:3.3f}$)".format(neighbors, eta), fontsize=10)

    ax.set_axis_off()
    ax.set_xlim(-0.5, 0.5)
    ax.set_ylim(0, 1)

    return


def surface_density_plot_no_yt(ax, snapnum, filename="keplerian_ring", density_limits=None, vlim=None):
    """
    Make the surface density plot (sans yt).

    Also returns the max and minimum values for the density so these can
    be passed to the next call, as well as vlim which are the colourmap
    max/min.
    """

    with h5py.File("{}_{:04d}.hdf5".format(filename, snapnum)) as filehandle:
        density = filehandle["PartType0"]["Density"][...]
        x, y = filehandle["PartType0"]["Coordinates"][:, 0:2].T

    new_vlim = (density.min(), density.max())

    if vlim is None:
        vlim = new_vlim

    ax.scatter(x, y, c=density, vmin=vlim[0], vmax=vlim[1], s=0.1)


    metadata = get_metadata("{}_{:04d}.hdf5".format(filename, snapnum))
    period = metadata["period"]

    t = ax.text(
        2.5,
        7.5,
        "Snapshot = {:04d}\nRotations = {:1.2f}".format(
            snapnum,
            float(metadata["header"]["Time"])/period
        ),
        color='black'
    )

    t.set_bbox(dict(alpha=0.5, color="white"))

    ax.axis("equal")

    ax.set_xlim(2, 8)
    ax.set_ylim(2, 8)

    # We now want to remove all of the ticklabels.

    for axis in ['x', 'y']:
        ax.tick_params(
            axis=axis,          
            which='both',      
            bottom='off',      
            top='off',         
            left='off',
            right='off',
            labelleft='off',
            labelbottom='off'
        ) 

    return density_limits, vlim



def surface_density_plot(ax, snapnum, filename="keplerian_ring", density_limits=None, vlim=None):
    """
    Make the surface density plot (via yt).

    Also returns the max and minimum values for the density so these can
    be passed to the next call, as well as vlim which are the colourmap
    max/min.
    """

    unit_base = {
        'length': (1.0, 'cm'),
        'velocity': (1.0, 'cm/s'),
        'mass': (1.0, 'g')
    }

    filename = "{}_{:04d}.hdf5".format(filename, snapnum)

    try:
        snap = yt.load(filename, unit_base=unit_base, over_refine_factor=2)
    except yt.utilities.exceptions.YTOutputNotIdentified:
        # Probably the file isn't here because we supplied a too high snapshot
        # number. Just return what we're given.
        return density_limits, vlim

    projection_plot = yt.ProjectionPlot(
        snap,
        "z",
        ("gas", "cell_mass"),
        width=5.5
    )

    max_density = snap.all_data()[("gas", "cell_mass")].max()
    min_density = snap.all_data()[("gas", "cell_mass")].min()
    
    new_density_limits = (min_density, max_density)

    if density_limits is None:
        density_limits = new_density_limits

    projection_plot.set_zlim("cell_mass", *density_limits)

    data = get_yt_actual_data(projection_plot, ("gas", "cell_mass"))

    # Becuase of the way plotting works, we also need a max/min for the colourmap.

    new_vlim = (data[0].min(), data[0].max())

    if vlim is None:
        vlim = new_vlim

    ax.imshow(
        data[0],
        cmap=data[1],
        vmin=vlim[0],
        vmax=vlim[1]
    )

    metadata = get_metadata(filename)
    period = metadata["period"]

    ax.text(
        20,
        80,
        "Snapshot = {:04d}\nRotations = {:1.2f}".format(
            snapnum,
            float(snap.current_time)/period
        ),
        color='white'
    )

    # We now want to remove all of the ticklabels.

    for axis in ['x', 'y']:
        ax.tick_params(
            axis=axis,          
            which='both',      
            bottom='off',      
            top='off',         
            left='off',
            right='off',
            labelleft='off',
            labelbottom='off'
        ) 

    return density_limits, vlim


def convert_snapshot_number_to_rotations_at(r, snapnum, filename):
    """
    Opens the file and extracts metadata to find the number of rotations.
    """

    metadata = get_metadata("{}_{:04d}.hdf5".format(filename, snapnum))

    t = metadata["period"]
    current_time = float(metadata["header"]["Time"])

    return current_time / t


if __name__ == "__main__":
    import argparse as ap

    parser = ap.ArgumentParser(
        description="""
                    Plotting code for the Keplerian Ring test. Uses yt to make
                    surface density plots.
                    """
    )

    parser.add_argument(
        "-b",
        "--beginning",
        help="""
             Initial snapshot to start analysis at.
             Default: First snapshot in the folder.
             """,
        default=-1,
        required=False
    )

    parser.add_argument(
        "-m",
        "--middle",
        help="""
             Middle snapshot in the top row of three surface density plots.
             Default: (end - beginning) // 2.
             """,
        default=-1,
        required=False
    )

    parser.add_argument(
        "-e",
        "--end",
        help="""
             Snapshot to end the analysis at.
             Default: Last snapshot in the folder.
             """,
        default=-1,
        required=False
    )

    parser.add_argument(
        "-f",
        "--filename",
        help="""
             Filename of the output.
             Default: plot.png
             """,
        default="plot.png",
        required=False
    )

    parser.add_argument(
        "-n",
        "--nbins",
        help="""
             Number of bins in histogramming (used for the density(r) plots).
             Default: 100.
             """,
        default=100,
        required=False
    )

    parser.add_argument(
        "-p",
        "--plotmassflow",
        help="""
             Plot the mass flow instead of several density(r) plots.
             Set this to a nonzero number to do this plot.
             Default: 0.
             """,
        default=0,
        required=False
    )

    parser.add_argument(
        "-s",
        "--slug",
        help="""
             The first part of the output filename. For example, for snapshots
             with the naming scheme keplerian_ring_xxxx.hdf5, this would be the
             default value of keplerian_ring.
             """,
        default="keplerian_ring",
        required=False
    )
    
    parser.add_argument(
        "-y",
        "--yt",
        help="""
             Use yt to do the plots at the top of the page. If set to anything
             other than a 'truthy' value, we will use a homebrew plotting
             setup. Default: False
             """,
        default=False,
        required=False
    )

    args = vars(parser.parse_args())

    filename = args["slug"]

    if args["beginning"] == -1:
        # We look for the maximum number of snapshots.
        numbers = [
            int(x[len(filename)+1:-5])
            for x in os.listdir() if
            x[:len(filename)] == filename and x[-1] == "5"
        ]

        snapshots = [min(numbers), (max(numbers) - min(numbers)) // 2, max(numbers)]
    else:
        snapshots = [args["beginning"], args["middle"], args["end"]]

    figure = plt.figure(figsize=(12, 10))
    axes = get_axes_grid(figure)

    density_limits = None
    vlim = None

    for snap, ax in zip(snapshots, tqdm(axes[0:3], desc="Images")):
        if args["yt"] and ytavail:
            density_limits, vlim = surface_density_plot(
                ax,
                snap,
                density_limits=density_limits,
                vlim=vlim
            )

            figure.subplots_adjust(hspace=0, wspace=0)
        else:
            density_limits, vlim = surface_density_plot_no_yt(
                ax,
                snap,
                density_limits=density_limits,
                vlim=vlim
            )

    # Now we need to do the density(r) plot.

    # Derived data includes density profiles and chi squared
    derived_data = get_derived_data(snapshots[0], snapshots[2], binnumber=int(args["nbins"]))

    if args["plotmassflow"]:
        plot_mass_flow(axes[3], snapshots[0], snapshots[2], derived_data[3])
    else:
        plot_density_r(axes[3], derived_data[0], derived_data[1], snapshots)

    plot_chisq(axes[4], snapshots[0], snapshots[2], derived_data[2])

    plot_extra_info(axes[5], "keplerian_ring_0000.hdf5")

    figure.savefig(args["filename"], dpi=300)


