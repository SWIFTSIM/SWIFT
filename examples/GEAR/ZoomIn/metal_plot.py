#!/usr/bin/env python3

import yt
import matplotlib.pyplot as plt


def save1DPlot(profiles):
    plt.figure(figsize=(8, 8))
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
    markers = ["s", "o"]
    for i, p in enumerate(profiles[0]):
        z = p.x.in_units("")
        m = p["Masses"].in_units("Msun").d
        plt.plot(z, m, linestyle="-", marker=markers[i],
                 markeredgecolor='none', linewidth=1.2, alpha=0.8)
    plt.semilogx()
    plt.semilogy()

    plt.xlabel(r"Metal fraction", fontsize='large')
    plt.ylabel("$\mathrm{Mass\ (Msun)}$", fontsize='large')
    plt.legend(profiles[1], loc=4, frameon=True, ncol=2, fancybox=True)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize='small')
    plt.grid(True)

    plt.savefig("metals.png", bbox_inches='tight', pad_inches=0.03, dpi=300)


def do1DPlot(f, name, i):
    sp = f.sphere(f.center, f.width)
    # Because ParticleProfilePlot doesn't exist, I will do the following trick.
    p = yt.create_profile(sp, ("PartType0", "Metallicity"),
                          ("PartType0", "Masses"),
                          weight_field=None, n_bins=50,
                          accumulation=False)

    return p
