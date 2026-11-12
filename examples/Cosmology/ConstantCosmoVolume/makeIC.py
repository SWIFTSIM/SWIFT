################################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

import h5py
import unyt
import swiftsimio as sw
import numpy as np

# Parameters
T_i = 100.0  # Initial temperature of the gas (in K)
z_i = 100.0  # Initial redshift
gamma = 5.0 / 3.0  # Gas adiabatic index
numPart_1D = 32
glassFile = "gravity_glassCube_32.hdf5"
fileName = "constantBox.hdf5"


# Some units
Mpc_in_m = 3.08567758e22
Msol_in_kg = 1.98848e30
Gyr_in_s = 3.08567758e19
mH_in_kg = 1.6737236e-27

# Some constants
kB_in_SI = 1.38064852e-23
G_in_SI = 6.67408e-11

# Some useful variables in h-full units
H_0 = 1.0 / Mpc_in_m * 10**5  # h s^-1
rho_0 = 3.0 * H_0**2 / (8 * np.pi * G_in_SI)  # h^2 kg m^-3
lambda_i = 64.0 / H_0 * 10**5  # h^-1 m (= 64 h^-1 Mpc)
x_min = -0.5 * lambda_i
x_max = 0.5 * lambda_i

# SI system of units
unit_l_in_si = Mpc_in_m
unit_m_in_si = Msol_in_kg * 1.0e10
unit_t_in_si = Gyr_in_s
unit_v_in_si = unit_l_in_si / unit_t_in_si
unit_u_in_si = unit_v_in_si**2

# ---------------------------------------------------

# Read the glass file
glass = sw.load(glassFile)

# Read particle positions and h from the glass
pos = glass.dark_matter.coordinates.value

# Total number of particles
numPart = np.size(pos) / 3
if numPart != numPart_1D**3:
    print("Non-matching glass file")
numPart = numPart_1D**3

# Set box size and interparticle distance
boxSize = x_max - x_min
delta_x = boxSize / numPart_1D

# Get the particle mass
a_i = 1.0 / (1.0 + z_i)
m_i = boxSize**3 * rho_0 / numPart

# Build the arrays
pos *= boxSize
v = np.zeros((numPart, 3))
ids = np.linspace(1, numPart, numPart)
m = np.zeros(numPart)
h = np.zeros(numPart)
u = np.zeros(numPart)

# Add in densities and material ids
rhos = np.ones(numPart) * rho_0
material_ids = np.zeros(numPart)

# Set the particles on the left
for i in range(numPart_1D):
    for j in range(numPart_1D):
        for k in range(numPart_1D):
            index = i * numPart_1D**2 + j * numPart_1D + k
            # coords[index,0] = (i + 0.5) * delta_x
            # coords[index,1] = (j + 0.5) * delta_x
            # coords[index,2] = (k + 0.5) * delta_x
            u[index] = kB_in_SI * T_i / (gamma - 1.0) / mH_in_kg
            h[index] = 1.2348 * delta_x
            m[index] = m_i
            v[index, 0] = 0.0
            v[index, 1] = 0.0
            v[index, 2] = 0.0

# Unit conversion
pos /= unit_l_in_si
v /= unit_v_in_si
m /= unit_m_in_si
h /= unit_l_in_si
u /= unit_u_in_si
rhos /= unit_m_in_si
rhos *= unit_l_in_si**3

boxSize /= unit_l_in_si

# File
UL = unyt.Unit(100.0 * unit_l_in_si * unyt.cm)
UM = unyt.Unit(1000.0 * unit_m_in_si * unyt.g)
Ut = unyt.Unit(1.0 * unit_t_in_si * unyt.s)
units = unyt.UnitSystem("constant_cosmo_volume_units", UL, UM, Ut)

w = sw.Writer(
    unit_system=units,
    boxsize=sw.cosmo_array(
        [boxSize, boxSize, boxSize],
        UL,
        comoving=True,
        scale_factor=1.0,
        scale_exponent=1,
    ),
    dimension=3,
)

w.gas.coordinates = pos * UL * w.a.comoving
w.gas.velocities = v * UL / Ut * w.a.comoving**0
w.gas.masses = m * UM * w.a.comoving**0
w.gas.smoothing_lengths = h * UL * w.a.comoving
w.gas.internal_energy = u * UL**2 / Ut**2 * w.a.comoving**0
w.gas.particle_ids = ids.astype(np.uint64)

w.write(fileName)

# Add the material IDs and densities
with h5py.File(fileName, "r+") as f:
    grp = f["PartType0"]
    grp.create_dataset("MaterialIDs", data=material_ids, dtype="i")
    grp.create_dataset("Density", data=rhos, dtype="f")
