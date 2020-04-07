###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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

import numpy as np
from astropy import units
from swiftsimio import Writer
import unyt

np.random.seed(50)

# parameters

filename = "simple_orbits.hdf5"
num_part = 1
masses = 1.
# If changed, need to update simple_orbits.yml
M = units.solMass.to("earthMass")
M = float(M)

min_r = 0.2
max_r = 5
boxsize = 2 * max_r

u_l = 1.49e13  # AU
u_m = 5.97e27  # Earth mass
u_v = 474047.  # AU / yr
u_t = u_l / u_v
G = 1.189972e-04  # grav. const.

# generate the coordinates
dist = np.random.rand(num_part) * (max_r - min_r) + min_r
angle = np.random.rand(num_part) * 2 * np.pi
# more easy to do it in 2D, therefore coords[:, 2] == 0
coords = np.zeros((num_part, 3))
coords[:, 0] = dist * np.sin(angle)
coords[:, 1] = dist * np.cos(angle)
coords += boxsize * 0.5

# generate the masses
m = np.ones(num_part) * masses

# generate the velocities
sign = np.random.rand(num_part)
sign[sign < 0.5] = -1
sign[sign >= 0.5] = 1

v = np.zeros((num_part, 3))
v[:, 0] = sign * np.sqrt(G * M / (dist * (1 + np.tan(angle)**2)))
v[:, 1] = - np.tan(angle) * v[:, 0]

# Write the snapshot
units = unyt.UnitSystem("Planets", unyt.AU, unyt.mearth, unyt.yr)

snapshot = Writer(units, boxsize * unyt.AU)
snapshot.dark_matter.coordinates = coords * unyt.AU
snapshot.dark_matter.velocities = v * unyt.AU / unyt.yr
snapshot.dark_matter.masses = m * unyt.mearth

snapshot.write(filename)
