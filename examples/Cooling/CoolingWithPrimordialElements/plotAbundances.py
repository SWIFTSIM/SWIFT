#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2023 Yves Revaz (yves.revaz@epfl.ch)
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

import matplotlib.pyplot as plt

plt.style.use("local.mplstyle")

import numpy as np
from glob import glob
import h5py
import sys

filenames = sys.argv[1:]

n = len(filenames)


rs = np.zeros(n)
nHIs = np.zeros(n)
nHIIs = np.zeros(n)
nHeIs = np.zeros(n)
nHeIIs = np.zeros(n)
nHeIIIs = np.zeros(n)
nHMs = np.zeros(n)
nH2Is = np.zeros(n)
nH2IIs = np.zeros(n)
nDIs = np.zeros(n)
nDIIs = np.zeros(n)
nHDIs = np.zeros(n)


for i, filename in enumerate(filenames):

    with h5py.File(filename, "r") as f:

        XHI = f["PartType0/HI"][:]
        XHII = f["PartType0/HII"][:]
        XHeI = f["PartType0/HeI"][:]
        XHeII = f["PartType0/HeII"][:]
        XHeIII = f["PartType0/HeIII"][:]
        XHM = f["PartType0/HM"][:]
        XH2I = f["PartType0/H2I"][:]
        XH2II = f["PartType0/H2II"][:]
        XDI = f["PartType0/DI"][:]
        XDII = f["PartType0/DII"][:]
        XHDI = f["PartType0/HDI"][:]

        redshift = f["Header"].attrs["Redshift"][0]
        rs[i] = redshift

        # compute hydrogen mass fraction
        XH = (
            XHI.mean()
            + XHII.mean()
            + XHM.mean()
            + XH2I.mean()
            + XH2II.mean()
            + XDI.mean()
            + XDII.mean()
            + XHDI.mean()
        )

        # make it equal to the hydrogen number density (do not consider density on purpose)
        nH = XH

        # all abundances per unit of nH

        nHI = XHI / nH
        nHII = XHII / nH
        nHeI = XHeI / nH / 4
        nHeII = XHeII / nH / 4
        nHeIII = XHeIII / nH / 4
        nHM = XHM / nH
        nH2I = XH2I / nH / 2
        nH2II = XH2II / nH / 2
        nDI = XDI / nH / 2
        nDII = XDII / nH / 2
        nHDI = XHDI / nH / 3

        nHIs[i] = nHI.mean()
        nHIIs[i] = nHII.mean()
        nHeIs[i] = nHeI.mean()
        nHeIIs[i] = nHeII.mean()
        nHeIIIs[i] = nHeIII.mean()
        nHMs[i] = nHM.mean()
        nH2Is[i] = nH2I.mean()
        nH2IIs[i] = nH2II.mean()
        nDIs[i] = nDI.mean()
        nDIIs[i] = nDII.mean()
        nHDIs[i] = nHDI.mean()


# data from Faure, A., Hily-Blant, Pineau des Forêts, Flower 2024
# The last lines have been removed on purpose. Their predictions
# are not reliable below z=20.
# Note: the "Redshift" field = z+1

with open("faure2024.dat", "r") as f:
    f.readline()
    f.readline()
    header = str.split(f.readline())
    lines = f.readlines()
    data = np.array(list(map(str.split, lines)), float)
    dic = {x: i for i, x in enumerate(header)}


##########################################
# plot
##########################################

plt.figure()
plt.subplots_adjust(left=0.17, bottom=0.12, right=0.97, top=0.97)

ax = plt.gca()

rs = rs + 1

pHI = ax.plot(rs, nHIs, ms=2, lw=1, label="HI")
pHII = ax.plot(rs, nHIIs, ms=2, lw=1, label="HII")
pHeI = ax.plot(rs, nHeIs, ms=2, lw=1, label="HeI")
pHeII = ax.plot(rs, nHeIIs, ms=2, lw=1, label="HeII")
pHeIII = ax.plot(rs, nHeIIIs, ms=2, lw=1, label="HeIII")
pHM = ax.plot(rs, nHMs, ms=2, lw=1, label="HM")
pH2I = ax.plot(rs, nH2Is, ms=2, lw=1, label="H2I")
pH2II = ax.plot(rs, nH2IIs, ms=2, lw=1, label="H2II")
pDI = ax.plot(rs, nDIs, ms=2, lw=1, label="DI")
pDII = ax.plot(rs, nDIIs, ms=2, lw=1, label="DII")
pHDI = ax.plot(rs, nHDIs, ms=2, lw=1, label="HDI")


# add data
ax.plot(data[:, dic["Redshift"]], data[:, dic["H"]], "--", color=pHI[0].get_color())
ax.plot(data[:, dic["Redshift"]], data[:, dic["H+"]], "--", color=pHII[0].get_color())
ax.plot(data[:, dic["Redshift"]], data[:, dic["He"]], "--", color=pHeI[0].get_color())
ax.plot(data[:, dic["Redshift"]], data[:, dic["He+"]], "--", color=pHeII[0].get_color())
ax.plot(data[:, dic["Redshift"]], data[:, dic["H-"]], "--", color=pHM[0].get_color())
ax.plot(data[:, dic["Redshift"]], data[:, dic["H2"]], "--", color=pH2I[0].get_color())
ax.plot(data[:, dic["Redshift"]], data[:, dic["H2+"]], "--", color=pH2II[0].get_color())
ax.plot(data[:, dic["Redshift"]], data[:, dic["D"]], "--", color=pDI[0].get_color())
ax.plot(data[:, dic["Redshift"]], data[:, dic["D+"]], "--", color=pDII[0].get_color())
ax.plot(data[:, dic["Redshift"]], data[:, dic["HD"]], "--", color=pHDI[0].get_color())


ax.set_xlabel(r"$\rm{1+z}$")
ax.set_xlim(100, 0.9)
ax.set_ylabel(r"$\rm{Fractional abundances}$")
ax.set_ylim(1e-21, 2)

ax.loglog()
ax.legend(loc="lower right")

plt.show()
