###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Stefan Arridge (stefan.arridge@durhama.ac.uk)
#                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
import unyt
import swiftsimio as sw

# Generates a SWIFT IC file with a constant density and pressure

# Parameters
boxSize = 1  # 1 kiloparsec
rho = 3.2e3  # Density in code units (3.2e6 is 0.1 hydrogen atoms per cm^3)
T = 4e3  # Initial Temperature
gamma = 5.0 / 3.0  # Gas adiabatic index
fileName = "coolingBox.hdf5"
# ---------------------------------------------------

# defines some constants
# need to be changed in plotTemperature.py too
h_frac = 0.76
mu = 4.0 / (1.0 + 3.0 * h_frac)

m_h_cgs = 1.67e-24
k_b_cgs = 1.38e-16

# defines units
unit_length = 3.0857e21  # kpc
unit_mass = 2.0e33  # solar mass
unit_time = 3.0857e16  # ~ Gyr

# Read id, position and h from glass
glass = sw.load("glassCube_32.hdf5")
ids = glass.gas.particle_ids.value
pos = glass.gas.coordinates.value * boxSize
h = glass.gas.smoothing_lengths.value * boxSize

# Compute basic properties
numPart = np.size(pos) // 3
mass = boxSize**3 * rho / numPart
internalEnergy = k_b_cgs * T * mu / ((gamma - 1.0) * m_h_cgs)
internalEnergy *= (unit_time / unit_length) ** 2

# File
UL = unyt.Unit(unit_length * unyt.cm)
UM = unyt.Unit(unit_mass * unyt.g)
Ut = unyt.Unit(unit_time * unyt.s)
units = unyt.UnitSystem("cooling_box_units", UL, UM, Ut)

w = sw.Writer(
    unit_system=units,
    boxsize=sw.cosmo_array(
        [boxSize] * 3, UL, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=3,
)

# Particle group
v = np.zeros((numPart, 3))
m = np.full(numPart, mass)
u = np.full(numPart, internalEnergy)

w.gas.coordinates = pos * UL * w.a.comoving
w.gas.velocities = v * UL / Ut * w.a.comoving**0
w.gas.masses = m * UM * w.a.comoving**0
w.gas.smoothing_lengths = h * UL * w.a.comoving
w.gas.internal_energy = u * UL**2 / Ut**2 * w.a.comoving**0
w.gas.particle_ids = ids.astype(np.uint64)

w.write(fileName)

print("Initial condition generated")
