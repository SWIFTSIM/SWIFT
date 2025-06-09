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


# ----------------------------------------------------
# Check that the total amount of ionized species
# remains constant
# ----------------------------------------------------

import os
import sys

import swiftsimio
from matplotlib import pyplot as plt
import argparse


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--basename", help="Snapshot basename", default="output")

    args = parser.parse_args()
    return args


def get_snapshot_list(snapshot_basename="output"):
    """
    Find the snapshot(s) that are to be plotted 
    and return their names as list
    """

    snaplist = []

    dirlist = os.listdir()
    for f in dirlist:
        if f.startswith(snapshot_basename) and f.endswith("hdf5"):
            snaplist.append(f)

    snaplist = sorted(snaplist)
    if len(snaplist) == 0:
        print(f"No snapshots with base '{snapshot_basename}' found")
        sys.exit(1)

    return snaplist


def compare_data(snaplist):
    """
    Create and save the plot
    """

    HI = []
    HII = []
    HeI = []
    HeII = []
    HeIII = []

    if "cosmo" in snaplist[0]:
        savename = "total_abundancies_cosmo.png"
    else:
        savename = "total_abundancies.png"

    for filename in snaplist:
        data = swiftsimio.load(filename)

        mXHI = data.gas.ion_mass_fractions.HI * data.gas.masses
        mXHII = data.gas.ion_mass_fractions.HII * data.gas.masses
        mXHeI = data.gas.ion_mass_fractions.HeI * data.gas.masses
        mXHeII = data.gas.ion_mass_fractions.HeII * data.gas.masses
        mXHeIII = data.gas.ion_mass_fractions.HeIII * data.gas.masses

        HI.append(mXHI.sum())
        HII.append(mXHII.sum())
        HeI.append(mXHeI.sum())
        HeII.append(mXHeII.sum())
        HeIII.append(mXHeIII.sum())

    plt.figure()
    plt.plot(range(len(snaplist)), HI, label="HI total mass")
    plt.plot(range(len(snaplist)), HII, label="HII total mass")
    plt.plot(range(len(snaplist)), HeI, label="HeI total mass")
    plt.plot(range(len(snaplist)), HeII, label="HeII total mass")
    plt.plot(range(len(snaplist)), HeIII, label="HeIII total mass")
    plt.legend()

    #  plt.show()
    plt.tight_layout()
    plt.savefig(savename, dpi=200)

    return


if __name__ == "__main__":
    # Get command line arguments
    args = parse_args()

    snapshot_base = args.basename
    snaplist = get_snapshot_list(snapshot_base)
    compare_data(snaplist)
