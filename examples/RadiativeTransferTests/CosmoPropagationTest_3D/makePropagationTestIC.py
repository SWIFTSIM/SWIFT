#!/usr/bin/env python3
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
#               2022 Tsang Keung Chan (chantsangkeung@gmail.com)
#               2024 Stan Verhoeve (s06verhoeve@gmail.com)
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
from swiftsimio import Writer
from swiftsimio.units import cosmo_units

glass = h5py.File("glassCube_64.hdf5", "r")
parts = glass["PartType0"]
xp = parts["Coordinates"][:]
h = parts["SmoothingLength"][:]
glass.close()

# replace the particle closest to the center
# by the star
r = np.sqrt(np.sum((0.5 - xp) ** 2, axis=1))
rmin = np.argmin(r)
mininds = np.argsort(r)
center_parts = xp[mininds[:4]]
xs = center_parts.sum(axis=0) / center_parts.shape[0]

# Double-check all particles for boundaries
for i in range(3):
    mask = xp[:, i] < 0.0
    xp[mask, i] += 1.0
    mask = xp[:, i] > 1.0
    xp[mask, i] -= 1.0

unitL = cosmo_units["length"]
edgelen = (2 * 260 * unyt.Mpc).to(unitL)
boxsize = unyt.unyt_array([edgelen.v, edgelen.v, edgelen.v], unitL)

xs = unyt.unyt_array(
    [np.array([xs[0] * edgelen, xs[1] * edgelen, xs[2] * edgelen])], unitL
)
xp *= edgelen
h *= edgelen

w = Writer(unit_system=cosmo_units, box_size=boxsize, dimension=3)

w.gas.coordinates = xp
w.stars.coordinates = xs
w.gas.velocities = np.zeros(xp.shape) * (unyt.cm / unyt.s)
w.stars.velocities = np.zeros(xs.shape) * (unyt.cm / unyt.s)
w.gas.masses = np.ones(xp.shape[0], dtype=float) * 1e1 * unyt.Msun
w.stars.masses = np.ones(xs.shape[0], dtype=float) * 100.0 * unyt.Msun
w.gas.internal_energy = (
    np.ones(xp.shape[0], dtype=float) * (300.0 * unyt.kb * unyt.K) / unyt.g
)

w.gas.smoothing_length = h
w.stars.smoothing_length = w.gas.smoothing_length[:1]

w.write("propagationTest-3D.hdf5")
