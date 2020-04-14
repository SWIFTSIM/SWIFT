#!/usr/bin/env python3

import yt
from yt.units import kpc, g, cm, K
from matplotlib.offsetbox import AnchoredText
import numpy as np

cmap_density = "algae"
cmap_temperature = "magma"
cmap_mass = "inferno"
cmap_metals = "coolwarm"
limits_density = (1e-10 * g / cm**2, 1e2 * g / cm**2)
limits_temperature = (10 * K, 1e5 * K)
limits_metals = None
limits_mass = None


def savePlot(fig, field):
    """
    Save the plot into a file.

    Parameters
    ----------

    fig: matplotlib figure
        The figure to save

    field: str
        The name of the field in the plot
    """
    fig.savefig("projection_%s.png" % field, bbox_inches="tight",
                pad_inches=0.03, dpi=300)


def doDensityPlot(f, name, i, fig, axes):
    """
    Generate a density projection plot.

    Parameters
    ----------

    f: yt dataset
        The data set to use in the projection

    name: str
        The name to print on the plot

    i: int
        The index of the plot

    fig: matplotlib figure
        The figure to use

    axes: list
        The list of axes to use for the plot
    """
    direction = "x"
    field = ("PartType0", "density")

    # compute the projection
    p = yt.ProjectionPlot(f, direction, field, center=f.center,
                          width=f.width, buff_size=(800, 800))

    # Compute the limits
    p.set_unit("density", "g / cm**2")
    data = p.to_fits_data()["density"].data
    if limits_density[0] > data.min():
        print("WARNING: data below min", data.min())
    if limits_density[1] < data.max():
        print("WARNING: data above max", data.max())

    # Adjust the plot
    p.set_cmap(field=field, cmap=cmap_density)
    p.set_log(field, True)
    p.set_zlim(field, limits_density[0], limits_density[1])

    # Draw it into the correct figure
    plot = p.plots[field]

    plot.figure = fig
    plot.axes = axes[i].axes

    # Add code name
    at = AnchoredText(name, loc=2, prop=dict(size=6), frameon=True)
    plot.axes.add_artist(at)

    # plot color bar
    if i != 0:
        plot.cax = axes.cbar_axes[0]
        p.hide_axes()
    p._setup_plots()

    if i == 0:
        z = 1. / f.scale_factor - 1.
        text = "Redshift = %.2g" % z
        prop = {
            "color": "w"
        }
        at = AnchoredText(text, loc=3, prop=prop, frameon=False)
        plot.axes.add_artist(at)

    # Add code name
    at = AnchoredText(name, loc=2, prop=dict(size=6), frameon=True)
    plot.axes.add_artist(at)


def doTemperaturePlot(f, name, i, fig, axes):
    """
    Generate a temperature projection plot.

    Parameters
    ----------

    f: yt dataset
        The data set to use in the projection

    name: str
        The name to print on the plot

    i: int
        The index of the plot

    fig: matplotlib figure
        The figure to use

    axes: list
        The list of axes to use for the plot
    """
    direction = "x"
    field = ("PartType0", "Temperature")
    a = f.scale_factor

    # compute the projection
    p = yt.ProjectionPlot(f, direction, field, center=f.center,
                          weight_field="density",
                          width=f.width, buff_size=(800, 800))

    # Compute the limits
    p.set_unit("Temperature", "K")
    data = p.to_fits_data()["Temperature"].data
    if limits_temperature[0] > data.min():
        print("WARNING: data below min", data.min())
    if limits_temperature[1] < data.max():
        print("WARNING: data above max", data.max())

    # Adjust the plot
    p.set_cmap(field=field, cmap=cmap_temperature)
    p.set_log(field, True)
    p.set_zlim(field, limits_temperature[0], limits_temperature[1])

    # Draw it into the correct figure
    plot = p.plots[field]

    plot.figure = fig
    plot.axes = axes[i].axes

    # Add code name
    at = AnchoredText(name, loc=2, prop=dict(size=6), frameon=True)
    plot.axes.add_artist(at)

    # plot color bar
    if i != 0:
        plot.cax = axes.cbar_axes[0]
        p.hide_axes()
    p._setup_plots()

    if i == 0:
        z = 1. / f.scale_factor - 1.
        text = "Redshift = %.2g" % z
        prop = {
            "color": "w"
        }
        at = AnchoredText(text, loc=3, prop=prop, frameon=False)
        plot.axes.add_artist(at)

    # Add code name
    at = AnchoredText(name, loc=2, prop=dict(size=6), frameon=True)
    plot.axes.add_artist(at)
