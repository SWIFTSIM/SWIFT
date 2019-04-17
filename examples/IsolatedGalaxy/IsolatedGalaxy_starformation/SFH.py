#!/usr/bin/env python3
################################################################################
# This file is part of SWIFT.
# Copyright (c) 2019 Folkert Nobels (nobels@strw.leidenuniv.nl)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################
import matplotlib

# matplotlib.use("Agg")
import h5py as h5
import numpy as np
import matplotlib.pyplot as plt

# check if we have tqdm installed
try:
    from tqdm import tqdm
except ImportError:
    raised_info = False

    def tqdm(x, *args, **kwargs):
        global raised_info

        if not raised_info:
            print("This script can display progress bars. Try `pip install tqdm`")
            raised_info = True
        return x


# Plot parameters
params = {
    "axes.labelsize": 20,
    "axes.titlesize": 20,
    "font.size": 16,
    "legend.fontsize": 16,
    "xtick.labelsize": 16,
    "ytick.labelsize": 16,
    "text.usetex": True,
    "figure.figsize": (7.5, 7.5),
    "figure.subplot.left": 0.12,
    "figure.subplot.right": 0.98,
    "figure.subplot.bottom": 0.08,
    "figure.subplot.top": 0.98,
    "lines.linewidth": 2.0,
    "text.latex.unicode": True,
}
plt.rcParams.update(params)
plt.rc("font", **{"family": "sans-serif", "sans-serif": ["Times"]})


def getSFH(filename, points):
    weightfac = 1e2 * points

    # Read the data
    with h5.File(filename, "r") as f:
        box_size = f["/Header"].attrs["BoxSize"][0]
        coordinates = f["/PartType4/Coordinates"][:, :]
        mass = f["/PartType4/Masses"][:]
        birth_time = f["/PartType4/BirthTime"][:]

    absmaxz = 2  # kpc
    absmaxxy = 10  # kpc

    part_mask = (
        ((coordinates[:, 0] - box_size / 2.0) > -absmaxxy)
        & ((coordinates[:, 0] - box_size / 2.0) < absmaxxy)
        & ((coordinates[:, 1] - box_size / 2.0) > -absmaxxy)
        & ((coordinates[:, 1] - box_size / 2.0) < absmaxxy)
        & ((coordinates[:, 2] - box_size / 2.0) > -absmaxz)
        & ((coordinates[:, 2] - box_size / 2.0) < absmaxz)
        & (birth_time > 0.0)
    )  # & (flag==1)

    birth_time = birth_time[part_mask]
    mass = mass[part_mask]

    histogram = np.histogram(
        birth_time, bins=points, weights=mass * weightfac, range=[0, 0.1]
    )
    values = histogram[0]
    xvalues = (histogram[1][:-1] + histogram[1][1:]) / 2.0
    return xvalues, values


def getsfrsnapwide(numbsnaps):
    """ Get the SFH from all the individual snaps
    """

    # time = np.arange(1, 101, 1)
    time = np.zeros(numbsnaps)
    SFR_gparticles = np.zeros(numbsnaps)
    previous_mass = 0
    previous_numb = 0
    for i in tqdm(range(0, numbsnaps)):
        # Read the data
        filename = "output_%04d.hdf5" % i
        with h5.File(filename, "r") as f:
            box_size = f["/Header"].attrs["BoxSize"][0]
            coordinates = f["/PartType0/Coordinates"][:, :]
            SFR = f["/PartType0/SFR"][:]
            coordinates_star = f["/PartType4/Coordinates"][:, :]
            masses_star = f["/PartType4/Masses"][:]
            time[i] = f["/Header"].attrs["Time"]

        absmaxz = 2  # kpc
        absmaxxy = 10  # kpc

        part_mask = SFR > 0

        SFR = SFR[part_mask]

        total_SFR = np.sum(SFR)
        SFR_gparticles[i] = total_SFR

    return time, SFR_gparticles


# run the main script
if __name__ == "__main__":
    # Read the logger file
    logdata = np.loadtxt("output_SFH_logger.txt")

    # Define the logger data in the correct units
    timelog = logdata[:, 1] * 9.778131e2

    # Calculate the cumulative sum of the elements of active sfh and formed stars
    CSFH_Mstar = np.cumsum(logdata[:, 4] * 1e10)
    CSFH_SFRdt = np.cumsum(logdata[:, 6] * 1e10)

    # plot the CSFH of the logger
    plt.plot(timelog, CSFH_Mstar, label="Stars formed")
    plt.plot(timelog, CSFH_SFRdt, label="Active gas particles")
    plt.xlabel("Time (Myr)")
    plt.ylabel("CSFH [$\\rm M_\odot$]")
    plt.xlim(0, 100)
    # plt.ylim(0, 1.2e9)
    plt.legend()
    plt.savefig("CSFH_logger.png")
    plt.close()

    # Calculate the Cumulative sum of the particles from the snap shots
    f = h5.File("./output_0100.hdf5", "r")
    birthtime = f["/PartType4/BirthTime"][:] * 9.778131e2
    mass = f["/PartType4/Masses"][:] * 1e10
    CSFH_birth = np.zeros(len(logdata[:, 0]))
    for i in tqdm(range(len(timelog))):
        mask = (birthtime > 0) & (birthtime <= timelog[i])
        CSFH_birth[i] = np.sum(mass[mask])

    # Plot the CSFH of the logger + from the birth time
    plt.plot(timelog, CSFH_Mstar, label="Stars formed")
    plt.plot(timelog, CSFH_SFRdt, label="Active gas particles")
    plt.plot(timelog, CSFH_birth, label="Birth time")
    plt.xlabel("Time (Myr)")
    plt.ylabel("CSFH [$\\rm M_\odot$]")
    plt.legend()
    plt.xlim(0, 100)
    plt.ylim(0, 1.2e9)
    plt.savefig("CSFH_all.png")
    plt.close()

    # Plot of the fractional difference between the different measures
    plt.yscale("log")
    plt.plot(
        timelog,
        np.abs(CSFH_Mstar - CSFH_SFRdt) / CSFH_Mstar,
        label="$\\frac{M_{log}-SFHdt_{log}}{M_{log}}$",
    )
    plt.plot(
        timelog,
        np.abs(CSFH_Mstar - CSFH_birth) / CSFH_Mstar,
        label="$\\frac{M_{log}-M_{birth}}{M_{birth}}$",
    )
    plt.plot(
        timelog,
        np.abs(CSFH_birth - CSFH_SFRdt) / CSFH_birth,
        label="$\\frac{M_{birth}-SFHdt_{log}}{M_{birth}}$",
    )
    plt.xlabel("Time (Myr)")
    plt.ylabel("Fractional difference")
    plt.xlim(0, 100)
    plt.ylim(1e-4, 1e0)
    plt.legend()
    plt.savefig("CSFH_fractional_diff.png")
    plt.close()

    # calculate the SFH from the last snapshot
    time_birth0, SFR_birth0 = getSFH("output_%04d.hdf5" % 100, 100)
    time_birth1, SFR_birth1 = getSFH("output_%04d.hdf5" % 100, 200)
    time_birth2, SFR_birth2 = getSFH("output_%04d.hdf5" % 100, 1000)
    time_birth3, SFR_birth3 = getSFH("output_%04d.hdf5" % 100, 4000)

    # make a plot of the different number of bins in the star formation routine
    plt.plot(time_birth3, SFR_birth3, label="4000 bins")
    plt.plot(time_birth2, SFR_birth2, label="1000 bins")
    plt.plot(time_birth1, SFR_birth1, label="200 bins")
    plt.plot(time_birth0, SFR_birth0, label="100 bins")
    plt.xlabel("Time (Myr)")
    plt.ylabel("SFH [$\\rm M_\odot \\rm yr^{-1}$]")
    plt.savefig("SFH_birth_time.png")
    plt.close()

    # Make a plot of the SFH from the snaps
    timesnap, SFRsnap = getsfrsnapwide(100)
    np.savetxt("SnapSFH.txt", np.transpose([timesnap, SFRsnap]))
    plt.plot(timesnap * 9.778131e2, SFRsnap * 1.022690e1, label="SFH gas tracers")
    plt.xlabel("Time (Myr)")
    plt.ylabel("SFH [$\\rm M_\odot \\rm yr^{-1}$]")
    plt.ylim(0, 15)
    plt.xlim(0, 100)
    plt.savefig("SFH_snapshots.png")
    plt.close()

    # Make a plot of the log file
    plt.plot(timelog, logdata[:, 7] * 1.023009e01, label="SFH log file")
    plt.xlabel("Time (Myr)")
    plt.ylabel("SFH [$\\rm M_\odot \\rm yr^{-1}$]")
    plt.ylim(0, 15)
    plt.xlim(0, 100)
    plt.legend()
    plt.savefig("SFH_log_file.png")
    plt.close()

    # Make a plot of the log file and the snaps
    plt.plot(timelog, logdata[:, 7] * 1.023009e01, label="SFH log file")
    plt.plot(timesnap * 9.778131e2, SFRsnap * 1.022690e1, label="SFH gas tracers")
    plt.xlabel("Time (Myr)")
    plt.ylabel("SFH [$\\rm M_\odot \\rm yr^{-1}$]")
    plt.ylim(0, 15)
    plt.xlim(0, 100)
    plt.legend()
    plt.savefig("SFH_all.png")
