#!/usr/bin/env python3

import yt
from yt.units import kpc
import matplotlib.pyplot as plt

width = 1500 * kpc


def save1DPlot(profiles):
    plt.figure(figsize=(8, 8))
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
    markers = ["s", "o"]
    for i, p in enumerate(profiles[0]):
        velocity = p.x.in_units("km/s").d
        mass = p["Masses"].in_units("Msun").d
        plt.plot(velocity, mass, linestyle="-", marker=markers[i],
                 markeredgecolor='none', linewidth=1.2, alpha=0.8)
    plt.semilogx()
    plt.semilogy()

    plt.xlabel("$\mathrm{Velocity\ (km/s)}$", fontsize='large')
    plt.ylabel(r"$\mathrm{Mass,}\/\mathrm{d}M\mathrm{/dlog}\/\mathrm{\rho}\/\mathrm{(M_{\odot})}$", fontsize='large')
    plt.legend(profiles[1], loc=4, frameon=True, ncol=2, fancybox=True)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize='small')
    plt.grid(True)

    plt.savefig("velocity.png", bbox_inches='tight', pad_inches=0.03, dpi=300)


def do1DPlot(f, name, i):
    sp = f.sphere("max", width)

    # Because ParticleProfilePlot doesn't exist, I will do the following trick.
    p = yt.create_profile(sp, ("PartType0", "velocity_magnitude"),  ("PartType0", "Masses"),
                          weight_field=None, n_bins=50,
                          accumulation=False)

    return p
