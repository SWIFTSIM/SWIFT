#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
from yt.units import Msun
from yt.utilities.cosmology import Cosmology
from multiprocessing import Pool
from functools import partial
from h5py import File

fe_sol_abund = 1.76603e-3
mg_sol_abund = 9.24316e-4


def saveSFRPlot(profiles):
    plt.figure(figsize=(8, 8))
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
    markers = ["s", "o"]

    sfr_fig = plt.figure()
    mstar_fig = plt.figure()

    # get the boundary
    t_max = None
    for p in profiles[0]:
        if t_max is None:
            t_max = p[1].max()
        else:
            tmp = p[1].max()
            if t_max < tmp:
                t_max = tmp
    t_max = float(t_max)

    for i, p in enumerate(profiles[0]):
        masses = p[0]
        form_time = p[1]
        time_range = [0, t_max]
        n_bins = 100
        hist, bins = np.histogram(form_time, bins=n_bins, range=time_range)
        inds = np.digitize(form_time, bins=bins)
        time = (bins[:-1] + bins[1:])/2

        sfr = np.array([masses[inds == j+1].sum()/(bins[j+1]-bins[j])
                        for j in range(len(time))])
        sfr[sfr == 0] = np.nan

        sfr /= 1e9  # convert M / Gyr -> M / yr

        plt.figure(sfr_fig.number)
        plt.plot(time, sfr, linestyle="-", marker=markers[i],
                 markeredgecolor="none", linewidth=1.2, alpha=0.8)

        mstars = np.array([masses[inds == j+1].sum()
                           for j in range(len(time))])
        mstars = np.cumsum(mstars) / 1e6
        plt.figure(mstar_fig.number)
        plt.plot(time, mstars, linestyle="-",
                 markeredgecolor="none", linewidth=1.2, alpha=0.8)

    plt.figure(sfr_fig.number)
    plt.xlabel('Time  [Gyr]')
    plt.ylabel('SFR  [M$_\odot$ yr$^{-1}$]')
    plt.legend(profiles[1])
    plt.savefig("sfr.png")

    plt.figure(mstar_fig.number)
    plt.xlabel('Time  [Gyr]')
    plt.ylabel('Stellar mass [$10^6$ M$_\odot$]')
    plt.legend(profiles[1])
    plt.savefig("mstars.png")


def doSFRPlot(f, name, i):
    sp = f.sphere(f.center, f.width)
    if name == "GEAR":
        masses = (sp["PartType1", "StarInitMass"] * 1e10 * Msun).in_units("Msun")
        a = sp["PartType1", "StarFormationTime"]
    else:
        masses = (sp["PartType4", "StarInitMass"] * 1e10 * Msun).in_units("Msun")
        a = sp["PartType4", "BirthScaleFactors"]
    cosmo = Cosmology(hubble_constant=f.hubble_constant,
                      omega_matter=f.omega_matter,
                      omega_lambda=f.omega_lambda)
    z = np.array(1. / a - 1.)
    with Pool() as p:
        formation_time = p.map(
            partial(cosmo.lookback_time, z_f=1e6), z)
    formation_time = f.arr(list(formation_time), "s").in_units("Gyr")

    return masses, formation_time


def saveAbundancesPlot(profiles):
    plt.figure(figsize=(8, 8))
    plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=0.5)
    markers = ["s", "o"]

    for i, p in enumerate(profiles[0]):
        mg = np.log10(p[0] / mg_sol_abund)
        fe = np.log10(p[1] / fe_sol_abund)

        plt.plot(fe, mg - fe, linestyle="", marker=markers[i],
                 markeredgecolor="none", alpha=0.8, markersize=3)

    plt.xlabel('[Fe / H]')
    plt.ylabel('[Mg / Fe]')
    plt.legend(profiles[1])
    plt.savefig("abundances.png")

    plt.xlim([-4, 0.5])
    plt.ylim([-1, 1.5])

    plt.savefig("abundances_zoom.png")


def getElement(f, sp, name, el):
    if name == "GEAR":
        metals = sp["PartType1", "StarMetals"]
        h = File(f.filename_template, "r")
        names = h["Header"].attrs["ChimieLabels"]
        names = names.split(b",")
        # remove garbage
        names = names[:-2]
        # bytes to unicode
        names = [name.decode() for name in names]
        h.close()
    else:
        metals = sp["PartType4", "ElementAbundances"]
        h = File(f.filename_template, "r")
        sub = h["SubgridScheme"].attrs
        names = []
        for i in range(metals.shape[1]):
            name = sub["Element %i" % i].decode()
            names.append(name)

    ind = names.index(el)
    return metals[:, ind]


def doAbundancesPlot(f, name, i):
    sp = f.sphere(f.center, f.width)

    mg = getElement(f, sp, name, "Mg")
    fe = getElement(f, sp, name, "Fe")

    return mg, fe
