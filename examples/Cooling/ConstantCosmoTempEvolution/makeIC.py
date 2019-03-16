################################################################################
# This file is part of SWIFT.
# Copyright (c) 2019 Stefan Arridge (stefan.arridge@durham.ac.uk)
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

from swiftsimio import Writer
from swiftsimio.units import cosmo_units

import unyt
import numpy as np
import h5py as h5


# Parameters
boxsize = 100 * unyt.Mpc
T_i = 1000.0  # Initial temperature of the gas (in K)
z_i = 100.0  # Initial redshift
gamma = 5.0 / 3.0  # Gas adiabatic index
n_p_1D = 32
glassFile = "gravity_glassCube_32.hdf5"
filename = "constantBox.hdf5"

# Cosmology (must be same as param file)
hubble_param = 0.6777  # same as in param file
Omega_bar = 0.0455  # same as in param file


# Read the glass file
glass = h5.File(glassFile, "r")

# Total number of particles
n_p = n_p_1D ** 3
# Read particle positions and from the glass
glass_pos = glass["/PartType1/Coordinates"][:, :]
glass.close()

# Calculate mean baryon density today from comological parameters
H_0 = 100.0 * hubble_param * unyt.km / unyt.s / unyt.Mpc
rho_crit_0 = 3.0 * H_0 ** 2 / (8.0 * np.pi * unyt.G)
rho_bar_0 = Omega_bar * rho_crit_0 

# From this, we calculate the mass of the gas particles
gas_particle_mass = rho_bar_0 * boxsize ** 3 / (n_p_1D ** 3)

# Generate object. cosmo_units corresponds to default Gadget-oid units
# of 10^10 Msun, Mpc, and km/s
x = Writer(cosmo_units, boxsize)

# 32^3 particles.
n_p = 32 ** 3

# Make gas coordinates from 0, 100 Mpc in each direction
x.gas.coordinates = glass_pos * boxsize

# Random velocities from 0 to 1 km/s
x.gas.velocities = np.zeros((n_p, 3)) * (unyt.km / unyt.s)

# Generate uniform masses as 10^6 solar masses for each particle
x.gas.masses = np.ones(n_p, dtype=float) * gas_particle_mass

# Generate internal energy corresponding to 10^4 K
x.gas.internal_energy = (
    np.ones(n_p, dtype=float) * (T_i * unyt.kb * unyt.K) / (1e6 * unyt.msun)
)

# Generate initial guess for smoothing lengths based on MIPS
x.gas.generate_smoothing_lengths(boxsize=boxsize, dimension=3)

# If IDs are not present, this automatically generates
x.write(filename)
