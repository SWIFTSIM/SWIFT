"""
Makes a rho-T plot. Uses the swiftsimio library.
"""

import matplotlib.pyplot as plt
import numpy as np

from swiftsimio import SWIFTDataset, SWIFTMetadata, SWIFTUnits

from unyt import mh, cm, Gyr
from tqdm import tqdm
from matplotlib.colors import LogNorm
from matplotlib.animation import FuncAnimation

# Constants; these could be put in the parameter file but are rarely changed.
density_bounds = [1e-8, 1e4]  # in nh/cm^3
temperature_bounds = [1e2, 1e8]  # in K
bins = 128

# Plotting controls
cmap = "viridis"


def get_data(filename):
    """
    Grabs the data (T in Kelvin and density in mh / cm^3).
    """

    data = SWIFTDataset(filename)

    data.gas.density.convert_to_units(mh / (cm ** 3))
    data.gas.temperature.convert_to_cgs()

    return data.gas.density, data.gas.temperature


def make_hist(filename, density_bounds, temperature_bounds, bins):
    """
    Makes the histogram for filename with bounds as lower, higher
    for the bins and "bins" the number of bins along each dimension.

    Also returns the edges for pcolormesh to use.
    """

    density_bins = np.logspace(
        np.log10(density_bounds[0]), np.log10(density_bounds[1]), bins
    )
    temperature_bins = np.logspace(
        np.log10(temperature_bounds[0]), np.log10(temperature_bounds[1]), bins
    )

    H, density_edges, temperature_edges = np.histogram2d(
        *get_data(filename), bins=[density_bins, temperature_bins]
    )

    return H.T, density_edges, temperature_edges


def setup_axes():
    """
    Creates the figure and axis object.
    """
    fig, ax = plt.subplots(1, figsize=(6, 5), dpi=300)

    ax.set_xlabel("Density [$n_H$ cm$^{-3}$]")
    ax.set_ylabel("Temperature [K]")

    ax.loglog()

    return fig, ax


def make_single_image(filename, density_bounds, temperature_bounds, bins):
    """
    Makes a single image and saves it to rhoTPlot_{filename}.png.
    
    Filename should be given _without_ hdf5 extension.
    """

    fig, ax = setup_axes()
    hist, d, T = make_hist(
        "{:s}.hdf5".format(filename), density_bounds, temperature_bounds, bins
    )

    mappable = ax.pcolormesh(d, T, hist, cmap=cmap, norm=LogNorm())
    fig.colorbar(mappable, label="Number of particles", pad=0)

    fig.tight_layout()

    fig.savefig("rhoTPlot_{:s}.png".format(filename))

    return


def make_movie(args, density_bounds, temperature_bounds, bins):
    """
    Makes a movie and saves it to rhoTPlot_{stub}.mp4.
    """

    fig, ax = setup_axes()

    def grab_metadata(n):
        filename = "{:s}_{:04d}.hdf5".format(args["stub"], n)
        data = SWIFTMetadata(filename)

        return data

    def grab_data(n):
        filename = "{:s}_{:04d}.hdf5".format(args["stub"], n)

        H, _, _ = make_hist(filename, density_bounds, temperature_bounds, bins)

        # Need to ravel because pcolormesh's set_array takes a 1D array. Might
        # as well do it here, beacuse 1d arrays are easier to max() than 2d.
        return H.ravel()

    histograms = [
        grab_data(n)
        for n in tqdm(
            range(args["initial"], args["final"] + 1), desc="Histogramming data"
        )
    ]

    metadata = [
        grab_metadata(n)
        for n in tqdm(
            range(args["initial"], args["final"] + 1), desc="Grabbing metadata"
        )
    ]

    units = SWIFTUnits("{:s}_{:04d}.hdf5".format(args["stub"], args["initial"]))

    # Need to get a reasonable norm so that we don't overshoot.
    max_particles = max([x.max() for x in histograms])

    norm = LogNorm(vmin=1, vmax=max_particles)

    # First, let's make the initial frame (we need this for our d, T values that we
    # got rid of in grab_data.
    hist, d, T = make_hist(
        "{:s}_{:04d}.hdf5".format(args["stub"], args["initial"]),
        density_bounds,
        temperature_bounds,
        bins,
    )

    mappable = ax.pcolormesh(d, T, hist, cmap=cmap, norm=norm)
    fig.colorbar(mappable, label="Number of particles", pad=0)

    fig.tight_layout()

    # Once we've rearranged the figure with tight_layout(), we can start laing
    # Down the metadata text.

    def format_metadata(metadata: SWIFTMetadata):
        t = metadata.t * units.units["Unit time in cgs (U_t)"]
        t.convert_to_units(Gyr)
        
        x = "$a$: {:2.2f}\n$z$: {:2.2f}\n$t$: {:2.2f}".format(
            metadata.a, metadata.z, t
        )

        return x

    text = ax.text(
        0.025,
        0.975,
        format_metadata(metadata[0]),
        ha="left",
        va="top",
        transform=ax.transAxes,
    )

    ax.text(
        0.975,
        0.975,
        metadata[0].code["Git Revision"].decode("utf-8"),
        ha="right",
        va="top",
        transform=ax.transAxes,
    )

    def animate(data):
        mappable.set_array(histograms[data])
        text.set_text(format_metadata(metadata[data]))

        return mappable

    animation = FuncAnimation(
        fig, animate, range(len(histograms)), fargs=[], interval=1000 / 25
    )

    animation.save("rhoTPlot_{:s}.mp4".format(args["stub"]))

    return


if __name__ == "__main__":
    import argparse as ap

    parser = ap.ArgumentParser(
        description="""
             Plotting script for making a rho-T plot.
             Takes the filename handle, start, and (optionally) stop
             snapshots. If stop is not given, png plot is produced for
             that snapshot. If given, a movie is made.
             """
    )

    parser.add_argument(
        "-i",
        "--initial",
        help="""Initial snapshot number. Default: 0.""",
        default=0,
        required=False,
        type=int,
    )

    parser.add_argument(
        "-f",
        "--final",
        help="""Final snapshot number. Default: 0.""",
        default=0,
        required=False,
        type=int,
    )

    parser.add_argument(
        "-s",
        "--stub",
        help="""Stub for the filename (e.g. santabarbara). This is
                the first part of the filename for the snapshots,
                not including the final underscore. Required.""",
        type=str,
        required=True,
    )

    args = vars(parser.parse_args())

    if args["final"] <= args["initial"]:
        # Run in single image mode.
        filename = "{:s}_{:04d}".format(args["stub"], args["initial"])

        make_single_image(
            filename,
            density_bounds=density_bounds,
            temperature_bounds=temperature_bounds,
            bins=bins,
        )

    else:
        # Movie mode!
        make_movie(
            args,
            density_bounds=density_bounds,
            temperature_bounds=temperature_bounds,
            bins=bins,
        )
