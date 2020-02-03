#!/usr/bin/env

from yt.extensions.astro_analysis.halo_analysis.api import HaloCatalog
import matplotlib.pyplot as plt
import numpy as np

Nbins = 20


def savePlot(data):
    plt.figure()
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
    markers = ["s", "o"]
    d_min = None
    d_max = None
    for d in data[0]:
        if d_min is None or d.min() < d_min:
            d_min = d.min()
        if d_max is None or d.max() > d_max:
            d_max = d.max()

    for i, d in enumerate(data[0]):
        plt.hist(d, Nbins, histtype="step", log=True, range=[d_min, d_max])

    plt.xlabel("$\mathrm{Mass\ \mathrm{(M_{\odot})}}$", fontsize='large')
    plt.ylabel(r"$\mathrm{Number}\ o\mathrm{f}\ \mathrm{Halos}$", fontsize='large')
    plt.legend(data[1], loc=4, frameon=True, ncol=2, fancybox=True)
    leg = plt.gca().get_legend()
    ltext = leg.get_texts()
    plt.setp(ltext, fontsize='small')
    plt.grid(True)

    plt.savefig("halo_distribution.png")


def doPlot(f, name, i):
    hc = HaloCatalog(data_ds=f, finder_method="fof")
    hc.create()

    masses = np.zeros(len(hc.catalog))
    for i, halo in enumerate(hc.catalog):
        masses[i] = halo["particle_mass"].in_units("Msun")

    return masses
