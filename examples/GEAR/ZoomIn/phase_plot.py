#!/usr/bin/env python3

import yt
from yt.units import Msun, amu, cm, K, kpc
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt
import matplotlib.lines as ln
import numpy as np

limits_temperature = (1e1 * K, 1e5 * K)
limits_density = (1e-7 * amu / cm**3, 1e7 * amu / cm**3)
limits_mass = (1e3 * Msun, 1e7 * Msun)


def save2DPlot(fig):
    fig.savefig("phase_2d.png", bbox_inches="tight",
                pad_inches=0.03, dpi=300)


def save1DPlotDensity(profiles):
    plt.figure(figsize=(8, 8))
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
    markers = ["s", "o"]
    for i, p in enumerate(profiles[0]):
        rho = p.x.in_units("g/cm**3").d
        mass = p["Masses"].in_units("Msun").d
        mass[mass == 0] = np.nan
        plt.plot(rho, mass, linestyle="-", marker=markers[i],
                 markeredgecolor='none', linewidth=1.2, alpha=0.8)
    plt.semilogx()
    plt.semilogy()

    plt.xlabel("$\mathrm{Density\ (g/cm^3)}$", fontsize='large')
    plt.ylabel(r"$\mathrm{Mass}\/\mathrm{(M_{\odot})}$", fontsize='large')
    plt.legend(profiles[1], loc=4, frameon=True, ncol=2, fancybox=True)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize='small')
    plt.grid(True)

    plt.savefig("phase_1d_density.png", bbox_inches='tight',
                pad_inches=0.03, dpi=300)


def save1DPlotTemperature(profiles):
    plt.figure(figsize=(8, 8))
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
    markers = ["s", "o"]
    for i, p in enumerate(profiles[0]):
        rho = p.x.in_units("K").d
        mass = p["Masses"].in_units("Msun").d
        mass[mass == 0] = np.nan
        plt.plot(rho, mass, linestyle="-", marker=markers[i],
                 markeredgecolor='none', linewidth=1.2, alpha=0.8)
    plt.semilogx()
    plt.semilogy()

    plt.xlabel("$\mathrm{Temperature\ K}$", fontsize='large')
    plt.ylabel(r"$\mathrm{Mass}\/\mathrm{(M_{\odot})}$", fontsize='large')
    plt.legend(profiles[1], loc=4, frameon=True, ncol=2, fancybox=True)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize='small')
    plt.grid(True)

    plt.savefig("phase_1d_temperature.png", bbox_inches='tight',
                pad_inches=0.03, dpi=300)


def do2DPlot(f, name, i, fig, axes):
    p = yt.PhasePlot(
        f, ("PartType0", "density"), ("PartType0", "Temperature"),
        ("PartType0", "mass"),
        weight_field=None, x_bins=300, y_bins=300)

    plt.title("Scale Factor %g" % f.scale_factor)

    # Set the unit system
    p.set_log("density", True)
    p.set_unit("density", "amu / cm**3")
    p.set_xlim(limits_density[0],
               limits_density[1])

    p.set_log("Temperature", True)
    p.set_unit("Temperature", "K")
    p.set_ylim(limits_temperature[0],
               limits_temperature[1])

    p.set_log("mass", True)
    p.set_unit("mass", "Msun")
    p.set_zlim("mass", limits_mass[0],
               limits_mass[1])

    plot = p.plots[("PartType0", "mass")]

    plot.figure = fig
    plot.axes = axes[i].axes

    # plot color bar
    if i != 0:
        plot.cax = axes.cbar_axes[0]
        # p.hide_axes()
    p._setup_plots()

    at = AnchoredText(name, loc=2, prop=dict(size=6), frameon=True)
    plot.axes.add_artist(at)

    # Make a grid
    x_ticks = plot.axes.get_xticks()
    y_ticks = plot.axes.get_yticks()
    for x in x_ticks:
        N = len(y_ticks)
        line = ln.Line2D([x]*N, y_ticks, linestyle=":", linewidth=0.6,
                         color="k", alpha=0.7)
        plot.axes.add_line(line)
    for y in y_ticks:
        N = len(x_ticks)
        line = ln.Line2D(x_ticks, [y]*N, linestyle=":", linewidth=0.6,
                         color="k", alpha=0.7)
        plot.axes.add_line(line)


def do1DPlotDensity(f, name, i):
    sp = f.sphere(f.center, f.width)
    # Because ParticleProfilePlot doesn't exist, I will do the following trick.
    p = yt.create_profile(sp, ("PartType0", "density"),
                          ("PartType0", "Masses"), weight_field=None,
                          n_bins=40, accumulation=False)
    return p


def do1DPlotTemperature(f, name, i):
    sp = f.sphere(f.center, f.width)
    # Because ParticleProfilePlot doesn't exist, I will do the following trick.
    p = yt.create_profile(sp, ("PartType0", "Temperature"),
                          ("PartType0", "Masses"), weight_field=None,
                          n_bins=40, accumulation=False)

    return p
