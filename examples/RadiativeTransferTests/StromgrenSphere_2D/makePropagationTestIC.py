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
# Intended for the photon propagation test.
# ---------------------------------------------------------------------

import h5py
import numpy as np
import unyt
import swiftsimio as sw

glass = h5py.File("glassPlane_128.hdf5", "r")
parts = glass["PartType0"]
xp = parts["Coordinates"][:]
h = parts["SmoothingLength"][:]
glass.close()

# find particles closest to the center
# and select a couple of them to put the star in their middle
r = np.sqrt(np.sum((0.5 - xp) ** 2, axis=1))
mininds = np.argsort(r)
center_parts = xp[mininds[:4]]
xs = center_parts.sum(axis=0) / center_parts.shape[0]

# Double-check all particles for boundaries
for i in range(2):
    mask = xp[:, i] < 0.0
    xp[mask, i] += 1.0
    mask = xp[:, i] > 1.0
    xp[mask, i] -= 1.0


unitL = unyt.cm
t_end = 1e-3 * unyt.s
edgelen = unyt.c.to("cm/s") * t_end * 2.0
edgelen = edgelen.to(unitL)

xs = unyt.unyt_array(
    [np.array([xs[0] * edgelen, xs[1] * edgelen, 0.0 * edgelen])], unitL
)
xp *= edgelen
h *= edgelen

boxsize_cosmo = sw.cosmo_array(
    [edgelen.v, edgelen.v], unitL, comoving=True, scale_factor=1.0, scale_exponent=1
)
w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system, boxsize=boxsize_cosmo, dimension=2
)

w.gas.coordinates = sw.cosmo_array(
    xp.value, xp.units, comoving=True, scale_factor=1.0, scale_exponent=1
)
w.stars.coordinates = sw.cosmo_array(
    xs.value, xs.units, comoving=True, scale_factor=1.0, scale_exponent=1
)
w.gas.velocities = sw.cosmo_array(
    np.zeros(xp.shape),
    unyt.cm / unyt.s,
    comoving=True,
    scale_factor=1.0,
    scale_exponent=0,
)
w.stars.velocities = sw.cosmo_array(
    np.zeros(xs.shape),
    unyt.cm / unyt.s,
    comoving=True,
    scale_factor=1.0,
    scale_exponent=0,
)
w.gas.masses = sw.cosmo_array(
    np.ones(xp.shape[0], dtype=np.float64) * 100,
    unyt.g,
    comoving=True,
    scale_factor=1.0,
    scale_exponent=0,
)
w.stars.masses = sw.cosmo_array(
    np.ones(xs.shape[0], dtype=np.float64) * 100,
    unyt.g,
    comoving=True,
    scale_factor=1.0,
    scale_exponent=0,
)
u = (300.0 * unyt.kb * unyt.K) / unyt.g
w.gas.internal_energy = sw.cosmo_array(
    np.ones(xp.shape[0], dtype=np.float64) * u.value,
    u.units,
    comoving=True,
    scale_factor=1.0,
    scale_exponent=-2,
)
w.gas.smoothing_lengths = sw.cosmo_array(
    h.value, h.units, comoving=True, scale_factor=1.0, scale_exponent=1
)
w.stars.smoothing_lengths = sw.cosmo_array(
    h.value[:1], h.units, comoving=True, scale_factor=1.0, scale_exponent=1
)

w.write("propagationTest-2D.hdf5")
