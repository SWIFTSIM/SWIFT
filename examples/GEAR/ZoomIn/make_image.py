#!/usr/bin/env python3

import yt
import sys
import projection_plot
import velocity_plot
import phase_plot
import halo_distribution
import metal_plot
import add_fields
import star_plot
import profile_plot
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
from yt.units import kpc
import numpy as np

# Parameters

swift = "swift_final.hdf5"
gear = "h050_final.hdf5"

width = 50 * kpc

do_dmo = True
do_hydro = True
do_stars = True
do_feedback = True
do_plot = {
    # DMO
    "projection_mass": do_dmo,
    "velocity": do_dmo,
    "halo_distribution": True,  # very slow
    "profile": do_dmo,
    # hydro
    "projection_density": do_hydro,
    "projection_temperature": do_hydro,
    "phase_1d_density": do_hydro,
    "phase_1d_temperature": do_hydro,
    "phase_2d": do_hydro,
    "metals": do_hydro,
    # stars
    "SFR": do_stars,
    "projection_mass_stars": do_stars,
    # feedback
    "abundances": do_feedback,
    "projection_metals": do_feedback,
    "iron_distribution": do_feedback,
    "mass_distribution": do_feedback
}

# Generate the figures

figsize = (100, 100)
figures = {
    "projection_density": plt.figure(figsize=figsize),
    "projection_temperature": plt.figure(figsize=figsize),
    "projection_metals": plt.figure(figsize=figsize),
    "projection_mass": plt.figure(figsize=figsize),
    "projection_mass_stars": plt.figure(figsize=figsize),
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
    "projection_mass_stars": AxesGrid(
        figures["projection_mass_stars"], fig_size, subgrid,
        add_all=True, share_all=True, cbar_mode="single",
        cbar_size="2%", cbar_pad=0.02),
    "phase_2d": AxesGrid(
        figures["phase_2d"], fig_size, subgrid, axes_pad=0.05,
        add_all=True, share_all=True, cbar_mode="single",
        cbar_size="2%", cbar_pad=0.05, aspect=False),
    "projection_metals": AxesGrid(
        figures["projection_metals"], fig_size, subgrid,
        add_all=True, share_all=True, cbar_mode="single",
        cbar_size="2%", cbar_pad=0.02),
}


# Data
data = {
    "phase_1d_density": ([], []),
    "phase_1d_temperature": ([], []),
    "halo_distribution": ([], []),
    "velocity": ([], []),
    "metals": ([], []),
    "SFR": ([], []),
    "profile": ([], []),
    "abundances": ([], []),
    "iron_distribution": ([], []),
    "mass_distribution": ([], []),
}

names = []


def savePlot():
    if do_plot["projection_density"]:
        projection_plot.savePlot(figures["projection_density"],
                                 "density")
    if do_plot["projection_temperature"]:
        projection_plot.savePlot(figures["projection_temperature"],
                                 "temperature")

    if do_plot["projection_metals"]:
        projection_plot.savePlot(figures["projection_metals"],
                                 "metals")

    if do_plot["projection_mass"]:
        projection_plot.savePlot(figures["projection_mass"],
                                 "mass")

    if do_plot["projection_mass_stars"]:
        projection_plot.savePlot(figures["projection_mass_stars"],
                                 "mass_stars")

    if do_plot["phase_1d_density"]:
        data["phase_1d_density"][1].extend(names)
        phase_plot.save1DPlotDensity(data["phase_1d_density"])

    if do_plot["phase_1d_temperature"]:
        data["phase_1d_temperature"][1].extend(names)
        phase_plot.save1DPlotTemperature(data["phase_1d_temperature"])

    if do_plot["phase_2d"]:
        phase_plot.save2DPlot(figures["phase_2d"])

    # halo distribution
    if do_plot["halo_distribution"]:
        data["halo_distribution"][1].extend(names)
        halo_distribution.savePlot(data["halo_distribution"])

    if do_plot["velocity"]:
        data["velocity"][1].extend(names)
        velocity_plot.save1DPlot(data["velocity"])

    if do_plot["metals"]:
        data["metals"][1].extend(names)
        metal_plot.save1DPlot(data["metals"])

    if do_plot["SFR"]:
        data["SFR"][1].extend(names)
        star_plot.saveSFRPlot(data["SFR"])

    if do_plot["profile"]:
        data["profile"][1].extend(names)
        profile_plot.savePlot(data["profile"])

    if do_plot["abundances"]:
        data["abundances"][1].extend(names)
        star_plot.saveAbundancesPlot(data["abundances"])

    if do_plot["iron_distribution"]:
        data["iron_distribution"][1].extend(names)
        star_plot.saveIronDistributionPlot(data["iron_distribution"])

    if do_plot["mass_distribution"]:
        data["mass_distribution"][1].extend(names)
        star_plot.saveMassDistributionPlot(data["mass_distribution"])


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

    if do_plot["projection_metals"]:
        projection_plot.doMetalsPlot(
            f, name, i, figures["projection_metals"],
            axes["projection_metals"])

    # Do mass projection plot
    if do_plot["projection_mass"]:
        projection_plot.doMassPlot(
            f, name, i, figures["projection_mass"],
            axes["projection_mass"], "all")

    # Do stellar mass projection plot
    if do_plot["projection_mass_stars"]:
        projection_plot.doMassPlot(
            f, name, i, figures["projection_mass_stars"],
            axes["projection_mass_stars"], "stars")

    # 1D Phase plot
    if do_plot["phase_1d_density"]:
        p = phase_plot.do1DPlotDensity(f, name, i)
        data["phase_1d_density"][0].append(p)

    if do_plot["phase_1d_temperature"]:
        p = phase_plot.do1DPlotTemperature(f, name, i)
        data["phase_1d_temperature"][0].append(p)

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

    # metal plot
    if do_plot["metals"]:
        p = metal_plot.do1DPlot(f, name, i)
        data["metals"][0].append(p)

    if do_plot["SFR"]:
        p = star_plot.doSFRPlot(f, name, i)
        data["SFR"][0].append(p)

    if do_plot["profile"]:
        p = profile_plot.doPlot(f, name, i)
        data["profile"][0].append(p)

    if do_plot["abundances"]:
        p = star_plot.doAbundancesPlot(f, name, i)
        data["abundances"][0].append(p)

    if do_plot["iron_distribution"]:
        p = star_plot.doIronDistributionPlot(f, name, i)
        data["iron_distribution"][0].append(p)

    if do_plot["mass_distribution"]:
        p = star_plot.doMassDistributionPlot(f, name, i)
        data["mass_distribution"][0].append(p)

    return center


center = None
# center = np.array([1724.33547783, 1802.56263082, 1785.09893269])
center = doPlot(gear, 1, "GEAR", center=center)
center = None
doPlot(swift, 0, "SWIFT", center=center)
savePlot()
