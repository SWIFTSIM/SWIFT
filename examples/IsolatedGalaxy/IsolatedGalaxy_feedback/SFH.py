"""
Plots the star formation history of the isolated galaxy, computed two ways:

 - from the summed star formation rates of the gas particles in each of the
   snapshots output_0001.hdf5 onwards;
 - from a mass-weighted histogram of the birth times of the star particles in
   the last snapshot.

Only particles within a 20x20x4 kpc slab centred on the box are used.
The two curves should agree if star formation behaves as expected.

Usage: python3 SFH.py
"""

import glob
import numpy as np
import matplotlib.pyplot as plt
import swiftsimio as sw

from tqdm import tqdm

# Half-sizes of the slab centred on the box, in kpc
absmaxxy = 10.0
absmaxz = 2.0


def in_slab(coordinates, boxsize):
    """
    Mask selecting the particles within the slab centred on the box.
    """

    offset = coordinates.to_value("kpc") - boxsize.to_value("kpc") / 2.0

    return (
        (np.abs(offset[:, 0]) < absmaxxy)
        & (np.abs(offset[:, 1]) < absmaxxy)
        & (np.abs(offset[:, 2]) < absmaxz)
    )


def getSFH(filename):
    """
    Star formation rate (in Msun/yr) from the birth times of the star
    particles, in 200 time bins (in Myr) from the start of the simulation to
    the time of the snapshot.
    """

    data = sw.load(filename)

    mask = in_slab(data.stars.coordinates, data.metadata.boxsize)
    birth_time = data.stars.birth_times[mask].to_value("Myr")
    mass = data.stars.masses[mask].to_value("Msun")

    bins = np.linspace(0.0, data.metadata.t.to_value("Myr"), 201)
    mass_formed, edges = np.histogram(birth_time, bins=bins, weights=mass)

    bin_width = (edges[1] - edges[0]) * 1e6  # yr
    xvalues = (edges[:-1] + edges[1:]) / 2.0
    return xvalues, mass_formed / bin_width


def getsfrsnapwide(filenames):
    """
    Time (in Myr) and summed star formation rate of the gas particles (in
    Msun/yr) of each snapshot.
    """

    time = np.zeros(len(filenames))
    SFR_gparticles = np.zeros(len(filenames))

    for i, filename in enumerate(tqdm(filenames)):
        data = sw.load(filename)

        SFR = data.gas.star_formation_rates.to_value("Msun/yr")
        part_mask = in_slab(data.gas.coordinates, data.metadata.boxsize) & (SFR > 0)

        time[i] = data.metadata.t.to_value("Myr")
        SFR_gparticles[i] = np.sum(SFR[part_mask])

    return time, SFR_gparticles


if __name__ == "__main__":

    snapshots = sorted(glob.glob("output_[0-9][0-9][0-9][0-9].hdf5"))

    time, SFR1 = getsfrsnapwide(snapshots[1:])
    time2, SFR3 = getSFH(snapshots[-1])
    plt.plot(time2[1:], SFR3[1:], label="Using birth_time of star particles")
    plt.plot(time, SFR1, label="Using SFR of gas particles", color="g")
    plt.xlabel("Time (Myr)")
    plt.ylabel(r"SFH ($\rm M_\odot \rm yr^{-1}$)")
    plt.ylim(0, 20)
    plt.legend()
    plt.savefig("SFH.png")
