#!/usr/bin/env python3

###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
##############################################################################


# ----------------------------------------------
# Plot the total photon energies in each photon
# group over the course of several snapshots
# ----------------------------------------------

import os
import sys

import matplotlib as mpl
import numpy as np
import swiftsimio
from matplotlib import pyplot as plt

# Parameters users should/may tweak
snapshot_base = "output"  # snapshot basename


# -----------------------------------------------------------------------

mpl.rcParams["text.usetex"] = True

# Read in cmdline arg: Are we plotting only one snapshot, or all?
plot_all = False  # plot all snapshots
try:
    snapnr = int(sys.argv[1])
except IndexError:
    plot_all = True


def get_snapshot_list(snapshot_basename="output"):
    """
    Find the snapshot(s) that are to be plotted 
    and return their names as list """

    snaplist = []

    if plot_all:
        dirlist = os.listdir()
        for f in dirlist:
            if f.startswith(snapshot_basename) and f.endswith("hdf5"):
                snaplist.append(f)

        snaplist = sorted(snaplist)

    else:
        fname = snapshot_basename + "_" + str(snapnr).zfill(4) + ".hdf5"
        if not os.path.exists(fname):
            print("Didn't find file", fname)
            quit(1)
        snaplist.append(fname)

    return snaplist


def plot_energies(snap_nrs, energy_sums):
    """
    Create the actual plot.
    """

    # Plot plot plot!
    fig = plt.figure(figsize=(5.0, 5.4), dpi=200)
    figname = "energy_budget.png"

    ax1 = fig.add_subplot(1, 1, 1)

    ngroups = energy_sums.shape[1]

    for g in range(ngroups):
        ax1.plot(snap_nrs, energy_sums[:, g], label=None)
        ax1.scatter(snap_nrs, energy_sums[:, g], label="group {0:d}".format(g + 1))

    ax1.set_xlabel("Snapshot")
    ax1.set_ylabel(
        r"Total energy [$" + energy_sums.units.latex_representation() + "$]",
        usetex=True,
    )

    # add title
    title = "Energy Budget"
    ax1.set_title(title)
    ax1.legend()

    plt.tight_layout()
    plt.savefig(figname)
    plt.close()

    return


def get_photon_energies(snaplist):
    """
    Get total photon energy in each photon group for a list of
    snapshots specified by `snaplist`

    snaplist: list of snapshot filenames

    returns:

        snap_nrs : list of integers of snapshot numbers
        energy_sums: np.array(shape=(len snaplist, ngroups)) of 
            total photon energies per group per snapshot
    """

    data = swiftsimio.load(snaplist[0])
    meta = data.metadata
    ngroups = int(meta.subgrid_scheme["PhotonGroupNumber"])

    energy_sums = np.zeros((len(snaplist), ngroups))
    snap_nrs = np.zeros(len(snaplist), dtype=int)

    for f, filename in enumerate(snaplist):

        data = swiftsimio.load(filename)

        for g in range(ngroups):
            en = getattr(data.gas.photon_energies, "group" + str(g + 1))
            energy_sums[f, g] = en.sum()

        nrstring = filename[len(snapshot_base) + 1 : -len(".hdf5")]
        nr = int(nrstring)
        snap_nrs[f] = nr

    energy_sums = energy_sums * en.units

    sortind = np.argsort(snap_nrs)
    snap_nrs = snap_nrs[sortind]
    energy_sums = energy_sums[sortind]

    return snap_nrs, energy_sums


if __name__ == "__main__":

    snaplist = get_snapshot_list(snapshot_base)
    snap_nrs, energy_sums = get_photon_energies(snaplist)

    plot_energies(snap_nrs, energy_sums)
