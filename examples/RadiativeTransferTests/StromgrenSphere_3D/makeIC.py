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

# ---------------------------------------------------------------------
# Add a single star in the center of a glass distribution
# The gas is set up with pure hydrogen gas.
# ---------------------------------------------------------------------

import h5py
import numpy as np
import unyt
import swiftsimio as sw
from swiftsimio.metadata.writer.unit_systems import cosmo_units

import stromgren_plotting_tools as spt

gamma = 5.0 / 3.0

# switch to replace the central gas particle with a star
# else put the star particle among gas particles
replace_gas = False


if __name__ == "__main__":

    glass = h5py.File("glassCube_64.hdf5", "r")
    parts = glass["PartType0"]
    xp = parts["Coordinates"][:]
    h = parts["SmoothingLength"][:]
    glass.close()

    r = np.sqrt(np.sum((0.5 - xp) ** 2, axis=1))

    if replace_gas:
        # replace a central gas particle with a star particle
        rmin = np.argmin(r)
        xs = xp[rmin]
        xp = np.delete(xp, rmin, axis=0)
        h = np.delete(h, rmin)
    else:
        # find particles closest to the center
        # and select a couple of them to put the star in their middle
        mininds = np.argsort(r)
        center_parts = xp[mininds[:4]]
        xs = center_parts.sum(axis=0) / center_parts.shape[0]

    # Double-check all particles for boundaries
    for i in range(3):
        mask = xp[:, i] < 0.0
        xp[mask, i] += 1.0
        mask = xp[:, i] > 1.0
        xp[mask, i] -= 1.0

    # Set up metadata
    unitL = unyt.Mpc
    edgelen = 22 * 1e-3 * unitL  # 22 so we can cut off 1kpc on each edge for image
    edgelen = edgelen.to(unitL)

    xs = unyt.unyt_array(
        [np.array([xs[0] * edgelen, xs[1] * edgelen, xs[2] * edgelen])], unitL
    )
    xp *= edgelen
    h *= edgelen

    boxsize_cosmo = sw.cosmo_array(
        [edgelen.value, edgelen.value, edgelen.value],
        edgelen.units,
        comoving=True,
        scale_factor=1.0,
        scale_exponent=1,
    )
    w = sw.Writer(unit_system=cosmo_units, boxsize=boxsize_cosmo, dimension=3)

    # write particle positions and smoothing lengths
    w.gas.coordinates = sw.cosmo_array(
        xp.value, xp.units, comoving=True, scale_factor=1.0, scale_exponent=1
    )
    w.stars.coordinates = sw.cosmo_array(
        xs.value, xs.units, comoving=True, scale_factor=1.0, scale_exponent=1
    )
    w.gas.velocities = sw.cosmo_array(
        np.zeros(xp.shape),
        unitL / unyt.Myr,
        comoving=True,
        scale_factor=1.0,
        scale_exponent=0,
    )
    w.stars.velocities = sw.cosmo_array(
        np.zeros(xs.shape),
        unitL / unyt.Myr,
        comoving=True,
        scale_factor=1.0,
        scale_exponent=0,
    )
    w.gas.smoothing_lengths = sw.cosmo_array(
        h.value, h.units, comoving=True, scale_factor=1.0, scale_exponent=1
    )
    w.stars.smoothing_lengths = sw.cosmo_array(
        h.value[:1], h.units, comoving=True, scale_factor=1.0, scale_exponent=1
    )

    # get gas masses
    XH = 1.0  # hydrogen mass fraction
    XHe = 0.0  # helium mass fraction
    nH = 1e-3 * unyt.cm ** (-3)
    rho_gas = nH * unyt.proton_mass / XH
    Mtot = rho_gas * edgelen**3
    mpart = Mtot / xp.shape[0]
    mpart = mpart.to(cosmo_units["mass"])
    w.gas.masses = sw.cosmo_array(
        np.ones(xp.shape[0], dtype=np.float64) * mpart.value,
        mpart.units,
        comoving=True,
        scale_factor=1.0,
        scale_exponent=0,
    )
    w.stars.masses = sw.cosmo_array(
        np.ones(xs.shape[0], dtype=np.float64) * mpart.value,
        mpart.units,
        comoving=True,
        scale_factor=1.0,
        scale_exponent=0,
    )

    # get gas internal energy for a given temperature and composition
    T = 100 * unyt.K
    XHI, XHII, XHeI, XHeII, XHeIII = spt.get_mass_fractions(T, XH, XHe)
    mu = spt.mean_molecular_weight(XHI, XHII, XHeI, XHeII, XHeIII)
    internal_energy = spt.internal_energy(T, mu, gamma)

    w.gas.internal_energy = sw.cosmo_array(
        np.ones(xp.shape[0], dtype=np.float64) * internal_energy.value,
        internal_energy.units,
        comoving=True,
        scale_factor=1.0,
        scale_exponent=-2,
    )

    w.write("stromgrenSphere-3D.hdf5")
