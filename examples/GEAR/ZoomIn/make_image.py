#!/usr/bin/env python3

import yt
from yt.units import Msun, kpc, s, km
import unyt
import sys
import matplotlib
import projection_plot
import phase_plot
import add_fields
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from matplotlib.offsetbox import AnchoredText
matplotlib.use('Agg')

# Parameters

snap = int(sys.argv[-1])

swift = "swift/snapshot_%04i.hdf5" % snap
gear = "gear/snapshot_%04i.hdf5" % snap


do_plot = {
    "projection_density": False,
    "projection_temperature": False,
    "projection_mass": True,
    "phase": False
}

# Generate the figures

figsize = (100, 100)
figures = {
    "projection_density": plt.figure(figsize=figsize),
    "projection_temperature": plt.figure(figsize=figsize),
    "projection_mass": plt.figure(figsize=figsize),
    "phase": plt.figure(figsize=figsize)
}

# define global variables
scale_factor = None

# Generate the axes

fig_size = (0.01, 0.01, 0.99, 0.99)
subgrid = (1, 2)
axes = {
    "projection_density": AxesGrid(
        figures["projection_density"], fig_size, subgrid,
        add_all=True, share_all=True, cbar_mode="single",
        cbar_size="2%", cbar_pad=0.02),
    "projection_temperature": AxesGrid(
        figures["projection_temperature"], fig_size, subgrid,
        add_all=True, share_all=True, cbar_mode="single",
        cbar_size="2%", cbar_pad=0.02),
    "projection_mass": AxesGrid(
        figures["projection_mass"], fig_size, subgrid,
        add_all=True, share_all=True, cbar_mode="single",
        cbar_size="2%", cbar_pad=0.02),
    "phase": AxesGrid(
        figures["phase"], fig_size, subgrid, axes_pad=0.05,
        add_all=True, share_all=True, cbar_mode="single",
        cbar_size="2%", cbar_pad=0.05, aspect=False)
}


def savePlot():
    if do_plot["projection_density"]:
        projection_plot.savePlot(figures["projection_density"],
                                 "density")
    if do_plot["projection_temperature"]:
        projection_plot.savePlot(figures["projection_temperature"],
                                 "temperature")

    if do_plot["projection_mass"]:
        projection_plot.savePlot(figures["projection_mass"],
                                 "mass")

    if do_plot["phase"]:
        phase_plot.savePlot(figures["phase"])


def doPlot(filename, i, name):
    f = yt.load(filename)
    if (do_plot["projection_temperature"] or do_plot["phase"]):
        add_fields.addTemperature(f)

    global scale_factor
    if scale_factor is None:
        scale_factor = f.scale_factor
    else:
        if abs(scale_factor - f.scale_factor) > scale_factor * 1e-2:
            raise Exception("Error: the snapshots are not at the same time")

    # Do density projection plot
    if do_plot["projection_density"]:
        projection_plot.doDensityPlot(
            f, name, i, figures["projection_density"],
            axes["projection_density"])

    # Do temperature projection plot
    if do_plot["projection_temperature"]:
        projection_plot.doTemperaturePlot(
            f, name, i, figures["projection_temperature"],
            axes["projection_temperature"])

    # Do mass projection plot
    if do_plot["projection_mass"]:
        add_fields.addMassDeposit(f)
        projection_plot.doMassPlot(
            f, name, i, figures["projection_mass"],
            axes["projection_mass"])

    # Phase plot
    if do_plot["phase"]:
        phase_plot.doPlot(f, name, i, figures["phase"],
                          axes["phase"])


doPlot(swift, 0, "SWIFT")
doPlot(gear, 1, "GEAR")
savePlot()
