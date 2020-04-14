#!/usr/bin/env python3

import yt
import projection_plot
import phase_plot
import add_fields
import star_plot
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from yt.units import kpc

swift = "swift_final.hdf5"
gear = "h050_final.hdf5"

width = 50 * kpc

do_dmo = True
do_hydro = True
do_stars = True
do_feedback = True
do_plot = {
    # hydro
    "projection_density": do_hydro,
    "projection_temperature": do_hydro,
    "phase_2d": do_hydro,
    # stars
    "SFR": do_stars,
    # feedback
    "abundances": do_feedback,
}

# Generate the figures

figsize = (100, 100)
figures = {
    "projection_density": plt.figure(figsize=figsize),
    "projection_temperature": plt.figure(figsize=figsize),
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
    "phase_2d": AxesGrid(
        figures["phase_2d"], fig_size, subgrid, axes_pad=0.05,
        add_all=True, share_all=True, cbar_mode="single",
        cbar_size="2%", cbar_pad=0.05, aspect=False),
}


# Data
data = {
    "SFR": ([], []),
    "abundances": ([], []),
}

names = []


def savePlot():
    if do_plot["projection_density"]:
        projection_plot.savePlot(figures["projection_density"],
                                 "density")
    if do_plot["projection_temperature"]:
        projection_plot.savePlot(figures["projection_temperature"],
                                 "temperature")

    if do_plot["phase_2d"]:
        phase_plot.save2DPlot(figures["phase_2d"])

    if do_plot["SFR"]:
        data["SFR"][1].extend(names)
        star_plot.saveSFRPlot(data["SFR"])

    if do_plot["abundances"]:
        data["abundances"][1].extend(names)
        star_plot.saveAbundancesPlot(data["abundances"])


def doPlot(filename, i, name, center):
    names.append(name)

    f = yt.load(filename)

    if center is None:
        center = f.find_max("Density")[1]

    f.center = center
    f.width = width * f.scale_factor

    if (do_plot["projection_temperature"] or do_plot["phase_2d"]):
        add_fields.addTemperature(f)
    add_fields.addMassDeposit(f)

    if (do_plot["projection_metals"] or do_plot["metals"]):
        add_fields.addMetals(f)

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

    # 2D Phase plot
    if do_plot["phase_2d"]:
        phase_plot.do2DPlot(f, name, i, figures["phase_2d"],
                            axes["phase_2d"])

    if do_plot["SFR"]:
        p = star_plot.doSFRPlot(f, name, i)
        data["SFR"][0].append(p)

    if do_plot["abundances"]:
        p = star_plot.doAbundancesPlot(f, name, i)
        data["abundances"][0].append(p)

    return center


center = None
# center = np.array([1724.33547783, 1802.56263082, 1785.09893269])
center = doPlot(gear, 1, "GEAR", center=center)
center = None
doPlot(swift, 0, "SWIFT", center=center)
savePlot()
