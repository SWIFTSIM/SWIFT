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
from getHMF import getHMFz, getHMFztinker

dlogm = 0.2
bins = 10 ** (np.arange(12, 15.2, dlogm))
V = 142.0 ** 3

itervalues = np.array([175, 185, 192, 198])

for j in itervalues:
    # Load the data
    g = h5py.File("stf_%04d.VELOCIraptor.properties" % j, "r")
    mass = g["Mass_200crit"][:] * 1e10  # convert to the correct unit
    binnedmass, massrange = np.histogram(mass, bins=bins)

    massnlarger = np.zeros(len(binnedmass))
    for i in range(0, len(massnlarger)):
        massnlarger[i] = np.sum(binnedmass[i:])

    f = h5py.File("snap_%04d.hdf5" % (j + 1))
    cosmo = f["Cosmology"]
    redshift = cosmo.attrs["Redshift"][0]
    a = cosmo.attrs["Scale-factor"][0]

    # Determine the HMF
    errormassn = massnlarger ** 0.5
    numbden = massnlarger / V / a ** 3
    numbdenerr = errormassn / V / a ** 3
    massplot = (massrange[0:15] + massrange[1:16]) / 2
    dernumbden = -np.diff(numbden) / np.diff(np.log10(massplot))
    dererr = 2 ** 0.5 / dlogm * (numbdenerr[0:14] + numbdenerr[1:15]) / 2

    plt.plot(
        (massplot[0:14] + massplot[1:15]) / 2, dernumbden, label="SWIFT - SPH $64^3$"
    )
    plt.fill_between(
        (massplot[0:14] + massplot[1:15]) / 2,
        dernumbden - dererr,
        dernumbden + dererr,
        alpha=0.4,
    )
    plt.xscale("log")
    plt.ylim(1e-6, 1e-1)
    plt.xlim(10 ** 11, 10 ** 15.5)

    xplace = 10 ** 14.5
    plt.text(xplace, 10 ** -2.3, "$\Omega_m=0.276$")
    plt.text(xplace, 10 ** -2.6, "$\Omega_b=0.0455$")
    plt.text(xplace, 10 ** -2.9, "$\Omega_\Lambda=0.724$")
    plt.text(xplace, 10 ** -3.2, "$h=0.703$")
    plt.text(xplace, 10 ** -3.5, "$z=%2.2f$" % redshift)

    m, dndlogm = getHMFz(redshift)
    plt.plot(m / 0.7, dndlogm * 0.7 ** 3, label="Sheth et al. 2001")

    m, dndlogm = getHMFztinker(redshift)
    plt.plot(m / 0.7, dndlogm * 0.7 ** 3, label="Tinker et al. 2008")

    plt.xlabel("M${}_{200}$ ($M_\odot$)")
    plt.ylabel("dn/d($\log$10(M${}_{200}$) ($Mpc^{-3}$)")
    plt.axvline(x=32 * 3.5e11, linestyle="--", color="k")
    plt.yscale("log")
    plt.legend()
    plt.savefig("./HMF_%04d.png" % j)
    plt.close()
