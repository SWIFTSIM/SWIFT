#!/usr/bin/env python3

import yt
from yt.units import Msun, amu, cm, K
from matplotlib.offsetbox import AnchoredText
import matplotlib.pyplot as plt


limits_temperature = (1e1 * K, 1e7 * K)
limits_density = (1e-7 * amu / cm**3, 1e7 * amu / cm**3)
limits_mass = (1e3 * Msun, 1e7 * Msun)

def savePlot(fig):
    fig.savefig("phase.png", bbox_inches="tight",
                pad_inches=0.03, dpi=300)


def doPlot(f, name, i, fig, axes):
    # limits in g / cm3
    # u_rho = g / cm**3
    # limits = [1e-33 * u_rho, 1e-24 * u_rho]

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
