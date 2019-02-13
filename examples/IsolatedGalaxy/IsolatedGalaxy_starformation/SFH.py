"""
Makes a movie using sphviewer and ffmpeg.

Edit low_frac and up_frac to focus on a certain view of the box.
The colour map can also be changed via colour_map.

Usage: python3 makeMovie.py CoolingHalo_

Written by James Willis (james.s.willis@durham.ac.uk)
"""

import glob
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

from tqdm import tqdm


def getSFH(filename):

    # Read the data
    with h5.File(filename, "r") as f:
        box_size = f["/Header"].attrs["BoxSize"][0]
        coordinates = f["/PartType4/Coordinates"][:, :]
        mass = f["/PartType4/Masses"][:]
        # flag = f["/PartType4/NewStarFlag"][:]
        birth_time = f["/PartType4/Birth_time"][:]

    absmaxz = 2  # kpc
    absmaxxy = 10  # kpc

    part_mask = (
        ((coordinates[:, 0] - box_size / 2.0) > -absmaxxy)
        & ((coordinates[:, 0] - box_size / 2.0) < absmaxxy)
        & ((coordinates[:, 1] - box_size / 2.0) > -absmaxxy)
        & ((coordinates[:, 1] - box_size / 2.0) < absmaxxy)
        & ((coordinates[:, 2] - box_size / 2.0) > -absmaxz)
        & ((coordinates[:, 2] - box_size / 2.0) < absmaxz)
    )  # & (flag==1)

    birth_time = birth_time[part_mask]
    mass = mass[part_mask]

    histogram = np.histogram(birth_time, bins=200, weights=mass * 2e4, range=[0, 0.1])
    values = histogram[0]
    xvalues = (histogram[1][:-1] + histogram[1][1:]) / 2.0
    return xvalues, values


def getsfrsnapwide():

    time = np.arange(1, 101, 1)
    SFR_sparticles = np.zeros(100)
    SFR_gparticles = np.zeros(100)
    new_sparticles = np.zeros(100)
    previous_mass = 0
    previous_numb = 0
    for i in tqdm(range(1, 100)):
        # Read the data
        filename = "output_%04d.hdf5" % i
        with h5.File(filename, "r") as f:
            box_size = f["/Header"].attrs["BoxSize"][0]
            coordinates = f["/PartType0/Coordinates"][:, :]
            SFR = f["/PartType0/SFR"][:]
            coordinates_star = f["/PartType4/Coordinates"][:, :]
            masses_star = f["/PartType4/Masses"][:]

        absmaxz = 2  # kpc
        absmaxxy = 10  # kpc

        part_mask = (
            ((coordinates[:, 0] - box_size / 2.0) > -absmaxxy)
            & ((coordinates[:, 0] - box_size / 2.0) < absmaxxy)
            & ((coordinates[:, 1] - box_size / 2.0) > -absmaxxy)
            & ((coordinates[:, 1] - box_size / 2.0) < absmaxxy)
            & ((coordinates[:, 2] - box_size / 2.0) > -absmaxz)
            & ((coordinates[:, 2] - box_size / 2.0) < absmaxz)
            & (SFR > 0)
        )

        SFR = SFR[part_mask]

        total_SFR = np.sum(SFR)
        SFR_gparticles[i] = total_SFR * 10

    return time[:-1], SFR_gparticles[1:]


if __name__ == "__main__":

    time, SFR1 = getsfrsnapwide()  # , SFR2, SFR_error = getsfrsnapwide()
    time2, SFR3 = getSFH("output_%04d.hdf5" % 100)
    plt.plot(time2[1:] * 1e3, SFR3[1:], label="Using birth_time of star particles")
    plt.plot(time, SFR1, label="Using SFR of gas particles", color="g")
    plt.xlabel("Time (Myr)")
    plt.ylabel("SFH ($\\rm M_\odot \\rm yr^{-1}$)")
    plt.ylim(0, 20)
    plt.legend()
    plt.savefig("SFH.png")
