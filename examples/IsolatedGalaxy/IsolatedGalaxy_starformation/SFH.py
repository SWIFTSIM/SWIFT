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
import glob
import numpy as np
import swiftsimio as sw
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


np.seterr(divide="ignore")
np.seterr(invalid="ignore")

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
}
plt.rcParams.update(params)


def getSFH(filename, points):
    """
    Star formation rate (in Msun/yr) from the birth times of the star
    particles within a 20x20x4 kpc slab centred on the box, in a number of
    time bins (in Myr) from the start of the simulation to the time of the
    snapshot.
    """

    data = sw.load(filename)

    offset = (
        data.stars.coordinates.to_value("kpc")
        - data.metadata.boxsize.to_value("kpc") / 2.0
    )
    birth_time = data.stars.birth_times.to_value("Myr")
    mass = data.stars.masses.to_value("Msun")

    absmaxz = 2  # kpc
    absmaxxy = 10  # kpc

    part_mask = (
        (np.abs(offset[:, 0]) < absmaxxy)
        & (np.abs(offset[:, 1]) < absmaxxy)
        & (np.abs(offset[:, 2]) < absmaxz)
        & (birth_time > 0.0)
    )

    birth_time = birth_time[part_mask]
    mass = mass[part_mask]

    bins = np.linspace(0.0, data.metadata.t.to_value("Myr"), points + 1)
    mass_formed, edges = np.histogram(birth_time, bins=bins, weights=mass)

    bin_width = (edges[1] - edges[0]) * 1e6  # yr
    xvalues = (edges[:-1] + edges[1:]) / 2.0
    return xvalues, mass_formed / bin_width


def getsfrsnapwide(filenames):
    """
    Get the SFH from all the individual snaps: time (in Myr) and summed star
    formation rate of the star-forming gas particles (in Msun/yr)
    """

    time = np.zeros(len(filenames))
    SFR_gparticles = np.zeros(len(filenames))
    for i, filename in enumerate(tqdm(filenames)):
        data = sw.load(filename)

        SFR = data.gas.star_formation_rates.to_value("Msun/yr")

        time[i] = data.metadata.t.to_value("Myr")
        SFR_gparticles[i] = np.sum(SFR[SFR > 0])

    return time, SFR_gparticles


# run the main script
if __name__ == "__main__":
    # Read the logger file
    logdata = np.loadtxt("SFR.txt")

    # The snapshots of the run
    snapshots = sorted(glob.glob("output_[0-9][0-9][0-9][0-9].hdf5"))

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
    data = sw.load(snapshots[-1])
    birthtime = data.stars.birth_times.to_value("Myr")
    mass = data.stars.masses.to_value("Msun")
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
    time_birth0, SFR_birth0 = getSFH(snapshots[-1], 100)
    time_birth1, SFR_birth1 = getSFH(snapshots[-1], 200)
    time_birth2, SFR_birth2 = getSFH(snapshots[-1], 1000)
    time_birth3, SFR_birth3 = getSFH(snapshots[-1], 4000)

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
    timesnap, SFRsnap = getsfrsnapwide(snapshots)
    np.savetxt(
        "SnapSFH.txt",
        np.transpose([timesnap, SFRsnap]),
        header="Time [Myr], SFR [Msun/yr]",
    )
    plt.plot(timesnap, SFRsnap, label="SFH gas tracers")
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
    plt.plot(timesnap, SFRsnap, label="SFH gas tracers")
    plt.xlabel("Time (Myr)")
    plt.ylabel("SFH [$\\rm M_\odot \\rm yr^{-1}$]")
    plt.ylim(0, 15)
    plt.xlim(0, 100)
    plt.legend()
    plt.savefig("SFH_all.png")
