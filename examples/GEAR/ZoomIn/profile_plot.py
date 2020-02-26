#!/usr/bin/env python3

import yt
from yt.units import Msun, kpc
import matplotlib.pyplot as plt

width = 20 * kpc
limits_mass = (1e3 * Msun, 1e7 * Msun)


def savePlot(profiles):
    plt.figure(figsize=(8, 8))
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
    markers = ["s", "o"]
    for i, p in enumerate(profiles[0]):
        r = p.x.in_units("pc").d
        mass = p["Masses"].in_units("Msun").d
        plt.plot(r, mass, linestyle="-", marker=markers[i],
                 markeredgecolor='none', linewidth=1.2, alpha=0.8)
    plt.semilogx()
    plt.semilogy()

    plt.xlabel(r"$\mathrm{Radius}\/\mathrm{(pc)}$", fontsize='large')
    plt.ylabel("$\mathrm{Integrated}\ \mathrm{mass}\ (Msun)}$",
               fontsize='large')
    plt.legend(profiles[1], loc=4, frameon=True, ncol=2, fancybox=True)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize='small')
    plt.grid(True)

    plt.savefig("integrated_mass.png", bbox_inches='tight',
                pad_inches=0.03, dpi=300)


def doPlot(f, name, i):
    sp = f.sphere(f.center, width)
    # Because ParticleProfilePlot doesn't exist, I will do the following trick.
    p = yt.create_profile(sp, "particle_radius", "Masses",
                          weight_field=None, n_bins=50,
                          accumulation=True)

    return p
