#!/usr/bin/env python3

import yt
from yt.units import Msun, kpc, s, km
import unyt
import sys
import matplotlib
import projection_plot
import velocity_plot
import phase_plot
import halo_distribution
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
    "projection_mass": False,
    "halo_distribution": False,
    "phase_1d": False,
    "phase_2d": False,
    "velocity": True
}

# Generate the figures

figsize = (100, 100)
figures = {
    "projection_density": plt.figure(figsize=figsize),
    "projection_temperature": plt.figure(figsize=figsize),
    "projection_mass": plt.figure(figsize=figsize),
    "phase_2d": plt.figure(figsize=figsize)
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
    "phase_2d": AxesGrid(
        figures["phase_2d"], fig_size, subgrid, axes_pad=0.05,
        add_all=True, share_all=True, cbar_mode="single",
        cbar_size="2%", cbar_pad=0.05, aspect=False)
}


# Data
data = {
    "phase_1d": ([], []),
    "halo_distribution": ([], []),
    "velocity": ([], []),
}

names = []


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

    if do_plot["phase_1d"]:
        data["phase_1d"][1].extend(names)
        phase_plot.save1DPlot(data["phase_1d"])

    if do_plot["phase_2d"]:
        phase_plot.save2DPlot(figures["phase_2d"])

    # halo distribution
    if do_plot["halo_distribution"]:
        data["halo_distribution"][1].extend(names)
        halo_distribution.savePlot(data["halo_distribution"])

    if do_plot["velocity"]:
        data["velocity"][1].extend(names)
        velocity_plot.save1DPlot(data["velocity"])


def doPlot(filename, i, name):
    names.append(name)

    f = yt.load(filename)
    if (do_plot["projection_temperature"] or do_plot["phase_2d"]):
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

    # 1D Phase plot
    if do_plot["phase_1d"]:
        p = phase_plot.do1DPlot(f, name, i)
        data["phase_1d"][0].append(p)

    # 2D Phase plot
    if do_plot["phase_2d"]:
        phase_plot.do2DPlot(f, name, i, figures["phase_2d"],
                            axes["phase_2d"])

    # halo distribution
    if do_plot["halo_distribution"]:
        m = halo_distribution.doPlot(f, name, i)
        data["halo_distribution"][0].append(m)

    # Velocity plot
    if do_plot["velocity"]:
        p = velocity_plot.do1DPlot(f, name, i)
        data["velocity"][0].append(p)


doPlot(swift, 0, "SWIFT")
doPlot(gear, 1, "GEAR")
savePlot()
