###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2023 Yves Revaz (yves.revaz@epfl.ch)
#                    Loic Hausammann (loic.hausammann@id.ethz.ch)
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

import h5py
import numpy as np
import unyt
import swiftsimio as sw

# Generates a SWIFT IC file with a constant density and pressure

# Parameters
periodic = 1  # 1 For periodic box
boxSize = 1  # 1 kiloparsec
# rho = 3.2e3  # Density in code units (3.2e6 is 0.1 hydrogen atoms per cm^3)

# Iliev test values
rho = 2471402.49842514  # Density in code units (1-hydrogen atoms per cm^3 assuming Hydrogen only)
T = 100  # Initial Temperature


gamma = 5.0 / 3.0  # Gas adiabatic index
fileName = "coolingBox.hdf5"
# ---------------------------------------------------

# defines some constants
# need to be changed in plotResults.py too
h_frac = 1.0  # 0.76
mu = 4.0 / (1.0 + 3.0 * h_frac)

m_h_cgs = 1.67262192369e-24
k_b_cgs = 1.380649e-16

# defines units
unit_length = 3.08567758e24  # Mpc
unit_mass = 1.98841e43  # 1e10 solar mass
unit_velocity = 1e5  # 1e5 km/s
unit_time = unit_length / unit_velocity

# Read id, position and h from glass
glass = h5py.File("glassCube_32.hdf5", "r")
ids = glass["/PartType0/ParticleIDs"][:]
pos = glass["/PartType0/Coordinates"][:, :] * boxSize
h = glass["/PartType0/SmoothingLength"][:] * boxSize

# Compute basic properties
numPart = np.size(pos) // 3
mass = boxSize**3 * rho / numPart
internalEnergy = k_b_cgs * T * mu / ((gamma - 1.0) * m_h_cgs)
internalEnergy *= (unit_time / unit_length) ** 2

# File
UL = unyt.Unit(unit_length * unyt.cm)
UM = unyt.Unit(unit_mass * unyt.g)
Ut = unyt.Unit(unit_time * unyt.s)
units = unyt.UnitSystem("cooling_heating_box_units", UL, UM, Ut)

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
