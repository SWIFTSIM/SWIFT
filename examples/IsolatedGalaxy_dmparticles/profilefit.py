#!/usr/bin/env python
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import BoundaryNorm
from matplotlib.ticker import MaxNLocator
import scipy.optimize as sco
import os


def linearfunc(x, a, b):
    return a * x + b


def radialfunc(r, h, A):
    return A * np.exp(-r / h) * r


def verticalfunc(z, A, z0, zoff):
    return 2 * A * np.exp(-(z - zoff) / z0)


def verticalfunc2(z, A, z0):
    return 2 * A * np.exp(-(z) / z0)


def verticalfunc3(z, A, z0, zoff, b):
    return 2 * A * np.exp(-(z - zoff) / z0) + b


Nmax = 2001
steps = 10
storefits = False
logfit = True
normalfit = False

# if the user wants to store the indivudal fits
if storefits:
    if not os.path.exists("radial"):
        os.mkdir("radial")
        os.mkdir("vertical")
        os.mkdir("histsnap")


# Initialize the arrays
R_ideal = np.linspace(0, 40, 100)
Z_ideal = np.linspace(0, 10, 100)

iterarray = np.arange(0, Nmax + 1, steps)

Z0t = np.zeros(len(iterarray))
Z0terr = np.zeros(len(iterarray))
h0t = np.zeros(len(iterarray))
h0terr = np.zeros(len(iterarray))
Ar = np.zeros(len(iterarray))
Arerr = np.zeros(len(iterarray))
Az = np.zeros(len(iterarray))
Azerr = np.zeros(len(iterarray))
time_array = np.zeros(len(iterarray))

ar = np.zeros(len(iterarray))
arerr = np.zeros(len(iterarray))
br = np.zeros(len(iterarray))
brerr = np.zeros(len(iterarray))
az = np.zeros(len(iterarray))
azerr = np.zeros(len(iterarray))
bz = np.zeros(len(iterarray))
bzerr = np.zeros(len(iterarray))
eps = 1e-6


for i in iterarray:
    # Getting the data from the snapshots
    f = h5py.File("output_%04d.hdf5" % i, "r")

    boxsize = f["Header"].attrs["BoxSize"] / 2.0

    time_array[int(i / steps)] = f["Header"].attrs["Time"]

    particles = f["PartType4"]
    coordinates = particles["Coordinates"][:, :]
    masses = particles["Masses"][:]

    R = (
        (coordinates[:, 0] - boxsize[0]) ** 2 + (coordinates[:, 1] - boxsize[1]) ** 2
    ) ** 0.5
    Z = np.abs(coordinates[:, 1] - boxsize[1])

    # Bin the coordinates to make them suitable for fitting
    Rhist = np.histogram(R, bins=100, range=[0, 40], normed=True)
    Zhist = np.histogram(Z, bins=100, range=[0, 10.0], normed=True)

    # Create correct variables for fitting
    Ry = Rhist[0]
    Rx = (Rhist[1][1:] + Rhist[1][: len(Rhist[0])]) / 2.0

    Zy = Zhist[0]
    Zx = (Zhist[1][1:] + Zhist[1][: len(Zhist[0])]) / 2.0

    # Fit with two methods: non-linear LSQ and linear LSQ in log space
    bestsolR = sco.curve_fit(radialfunc, Rx[10:], Ry[10:], p0=[2.0, 0.2])
    bestsolZ = sco.curve_fit(verticalfunc, Zx[40:], Zy[40:])
    bestsolRlog = sco.curve_fit(linearfunc, Rx[10:], np.log10(Ry[10:] + eps))
    bestsolZlog = sco.curve_fit(linearfunc, Zx[40:], np.log10(Zy[40:] + eps))

    # Store variables
    h0t[int(i / steps)] = bestsolR[0][0]
    Z0t[int(i / steps)] = bestsolZ[0][1]
    Ar[int(i / steps)] = bestsolR[0][1]
    Az[int(i / steps)] = bestsolZ[0][0]
    Z0terr[int(i / steps)] = (bestsolZ[1][1, 1]) ** 0.5
    h0terr[int(i / steps)] = (bestsolR[1][0, 0]) ** 0.5
    Arerr[int(i / steps)] = (bestsolR[1][1, 1]) ** 0.5
    Azerr[int(i / steps)] = (bestsolZ[1][0, 0]) ** 0.5

    ar[int(i / steps)] = bestsolRlog[0][0]
    arerr[int(i / steps)] = (bestsolRlog[1][0, 0]) ** 0.5
    br[int(i / steps)] = bestsolRlog[0][1]
    brerr[int(i / steps)] = (bestsolRlog[1][1, 1]) ** 0.5
    az[int(i / steps)] = bestsolZlog[0][0]
    azerr[int(i / steps)] = (bestsolZlog[1][0, 0]) ** 0.5
    bz[int(i / steps)] = bestsolZlog[0][1]
    bzerr[int(i / steps)] = (bestsolZlog[1][1, 1]) ** 0.5

    if storefits:
        plt.step(Rx, Ry)
        plt.plot(
            R_ideal,
            radialfunc(R_ideal, bestsolR[0][0], bestsolR[0][1]),
            label="Non linear LSQ",
        )
        plt.plot(
            R_ideal,
            10 ** (linearfunc(R_ideal, bestsolRlog[0][0], bestsolRlog[0][1])),
            label="Linear LSQ",
        )
        plt.xlim(0, 40)
        plt.ylim(0, 0.25)
        plt.xlabel("R (kpc)")
        plt.ylabel("Probability")
        plt.savefig("./radial/radialsnap%04d.png" % i)
        plt.close()

        plt.step(Zx, Zy)
        plt.plot(
            Z_ideal,
            verticalfunc(Z_ideal, bestsolZ[0][0], bestsolZ[0][1], bestsolZ[0][2]),
            label="Non linear LSQ",
        )
        plt.plot(
            Z_ideal,
            10 ** (linearfunc(Z_ideal, bestsolZlog[0][0], bestsolZlog[0][1])),
            label="Linear LSQ",
        )
        plt.xlim(0, 10.0)
        plt.ylim(0, 0.6)
        plt.xlabel("z (kpc)")
        plt.ylabel("Probability")
        plt.savefig("./vertical/verticalsnap%04d.png" % i)
        plt.close()

time_array[-1] = 2.0

ax = plt.subplot(111)
ax.set_yscale("log")
if logfit:
    plt.errorbar(
        time_array,
        np.absolute(az / (az[0]) - 1),
        yerr=azerr / (az[0]),
        label="z0 scale height (Log space)",
    )
    plt.errorbar(
        time_array,
        np.absolute(ar / (ar[0]) - 1),
        yerr=arerr / (ar[0]),
        label="h scale lenght (Log space)",
    )
if normalfit:
    plt.errorbar(
        time_array,
        np.absolute(Z0t / (Z0t[0]) - 1),
        yerr=Z0terr / (Z0t[0]),
        label="z0 scale height (normal space)",
    )
    plt.errorbar(
        time_array,
        np.absolute(h0t / (h0t[0]) - 1),
        yerr=h0terr / (h0t[0]),
        label="h scale height (normal space)",
    )
ax.set_xlabel("Time (Gyr)")
ax.set_ylabel("Fractional difference")
plt.legend()
plt.savefig("Fitdifference-witherror.pdf")
plt.close()


ax = plt.subplot(111)
ax.set_yscale("log")
if logfit:
    plt.plot(
        time_array, np.absolute(az / (az[0]) - 1), label="z0 scale height (Log space)"
    )
    plt.plot(
        time_array, np.absolute(ar / (ar[0]) - 1), label="h scale lenght (Log space)"
    )
if normalfit:
    plt.plot(
        time_array,
        np.absolute(Z0t / (Z0t[0]) - 1),
        label="z0 scale height (normal space)",
    )
    plt.plot(
        time_array,
        np.absolute(h0t / (h0t[0]) - 1),
        label="h scale height (normal space)",
    )
ax.set_xlabel("Time (Gyr)")
ax.set_ylabel("Fractional difference")
plt.legend()
plt.savefig("Fitdifference.pdf")
plt.show()
