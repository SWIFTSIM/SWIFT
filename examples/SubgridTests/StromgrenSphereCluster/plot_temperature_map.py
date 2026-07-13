"""
Makes a mass-weighted gas temperature projection plot. Reuses the
get_gas_temperatures()/get_gas_mu() helpers and the mass-projection plotting
conventions from plot_gas_density.py in this same directory (same cooling
model dispatch, same star/box-overlay style) -- the only new piece is the
temperature projection itself, since project_gas() gives a mass-weighted
COLUMN of whatever quantity you feed it, so a true (density-independent)
average temperature per pixel needs num = project(mass*T), den =
project(mass), then num/den.
"""

import matplotlib.pyplot as plt
import numpy as np

import swiftsimio as sw
from swiftsimio.visualisation.projection import project_gas

import unyt
from unyt import mh, cm
from matplotlib.colors import LogNorm

from plot_gas_density import get_gas_temperatures, get_sink_and_stars_positions

image_resolution = 1024
cmap = "inferno"


def make_temperature_projection(filename, image_resolution):
    data = sw.load(filename)

    data.gas.temperatures = get_gas_temperatures(data)
    data.gas.mass_weighted_temperature = data.gas.masses * data.gas.temperatures

    numerator = project_gas(
        data,
        resolution=image_resolution,
        project="mass_weighted_temperature",
        parallel=True,
        periodic=True,
    )
    denominator = project_gas(
        data,
        resolution=image_resolution,
        project="masses",
        parallel=True,
        periodic=True,
    )

    # Avoid 0/0 in empty pixels; mask them out instead of showing bogus T.
    with np.errstate(invalid="ignore", divide="ignore"):
        mean_temperature = (numerator / denominator).to(unyt.K)
    mean_temperature = np.where(
        denominator.value > 0, mean_temperature.value, np.nan
    )

    boxsize = data.metadata.boxsize
    x_edges = np.linspace(0 * unyt.kpc, boxsize[0], image_resolution)
    y_edges = np.linspace(0 * unyt.kpc, boxsize[1], image_resolution)

    return mean_temperature.T, x_edges, y_edges


def make_single_image(filename, image_resolution):
    file = "{:s}.hdf5".format(filename)
    fig, ax = plt.subplots(1, figsize=(6, 5), dpi=300)
    ax.set_xlabel("x [kpc]")
    ax.set_ylabel("y [kpc]")

    T_map, x_edges, y_edges = make_temperature_projection(file, image_resolution)

    mappable = ax.pcolormesh(
        x_edges, y_edges, T_map, cmap=cmap, norm=LogNorm(vmin=10, vmax=2e4)
    )
    fig.colorbar(mappable, label="Mass-weighted temperature [K]", pad=0)

    sink_pos, star_pos = get_sink_and_stars_positions(file)
    if star_pos.size != 0:
        ax.scatter(star_pos[:, 0], star_pos[:, 1], c="limegreen", zorder=1, marker="*")

    ax.text(
        0.7, 0.95, "$N_{\\mathrm{star}}" + " = {}$".format(len(star_pos)),
        transform=ax.transAxes, fontsize=7,
        bbox=dict(facecolor="white", alpha=0.8),
    )

    fig.tight_layout()
    image = "temperature_projection_{:s}.png".format(filename[-4:])
    fig.savefig(image)
    return


if __name__ == "__main__":
    import argparse as ap

    parser = ap.ArgumentParser(
        description="Plotting script for a mass-weighted temperature projection."
    )
    parser.add_argument("-i", "--initial", default=0, type=int)
    parser.add_argument("-s", "--stub", type=str, required=True)
    args = vars(parser.parse_args())

    filename = "{:s}_{:04d}".format(args["stub"], args["initial"])
    make_single_image(filename, image_resolution=image_resolution)
