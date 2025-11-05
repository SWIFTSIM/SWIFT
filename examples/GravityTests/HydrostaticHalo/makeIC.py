###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Stefan Arridge (stefan.arridge@durham.ac.uk)
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
import sys
import numpy as np
import math
import os

# Generates N particles in a spherically symmetric distribution with density profile ~r^(-2)
# usage: python3 makeIC.py 1000: generate 1000 particles

# Some constants
OMEGA = 0.3  # Cosmological matter fraction at z = 0
PARSEC_IN_CGS = 3.0856776e18
KM_PER_SEC_IN_CGS = 1.0e5
CONST_G_CGS = 6.672e-8
h = 0.67777  # hubble parameter
gamma = 5.0 / 3.0
eta = 1.2349

# First set unit velocity and then the circular velocity parameter for the isothermal potential
const_unit_velocity_in_cgs = 1.0e5  # kms^-1
v_c = 200.0
v_c_cgs = v_c * const_unit_velocity_in_cgs

# Now we use this to get the virial mass and virial radius, which we will set to be the unit mass and radius

# Find H_0, the inverse Hubble time, in cgs
H_0_cgs = 100.0 * h * KM_PER_SEC_IN_CGS / (1.0e6 * PARSEC_IN_CGS)

# From this we can find the virial radius, the radius within which the average density of the halo is
# 200. * the mean matter density
r_vir_cgs = v_c_cgs / (10.0 * H_0_cgs * np.sqrt(OMEGA))

# Now get the virial mass
M_vir_cgs = r_vir_cgs * v_c_cgs ** 2 / CONST_G_CGS

# Now set the unit length and mass
const_unit_mass_in_cgs = M_vir_cgs
const_unit_length_in_cgs = r_vir_cgs


# derived quantities
const_unit_time_in_cgs = const_unit_length_in_cgs / const_unit_velocity_in_cgs
const_G = (
    CONST_G_CGS
    * const_unit_mass_in_cgs
    * const_unit_time_in_cgs
    * const_unit_time_in_cgs
    / (const_unit_length_in_cgs * const_unit_length_in_cgs * const_unit_length_in_cgs)
)

# Parameters
periodic = 1  # 1 For periodic box
boxSize = 4.0
G = const_G
N = int(sys.argv[1])  # Number of particles

# set seed for random number
np.random.seed(1234)

# Positions
# r^(-2) distribution corresponds to uniform distribution in radius
radius = (
    boxSize * np.sqrt(3.0) / 2.0 * np.random.rand(N)
)  # the diagonal extent of the cube
ctheta = -1.0 + 2 * np.random.rand(N)
stheta = np.sqrt(1.0 - ctheta ** 2)
phi = 2 * math.pi * np.random.rand(N)
coords = np.zeros((N, 3))
coords[:, 0] = radius * stheta * np.cos(phi)
coords[:, 1] = radius * stheta * np.sin(phi)
coords[:, 2] = radius * ctheta

# shift to centre of box
coords += np.full((N, 3), boxSize / 2.0)

# now find the particles which are within the box

x_coords = coords[:, 0]
y_coords = coords[:, 1]
z_coords = coords[:, 2]

ind = np.where(x_coords < boxSize)[0]
x_coords = x_coords[ind]
y_coords = y_coords[ind]
z_coords = z_coords[ind]

ind = np.where(x_coords > 0.0)[0]
x_coords = x_coords[ind]
y_coords = y_coords[ind]
z_coords = z_coords[ind]

ind = np.where(y_coords < boxSize)[0]
x_coords = x_coords[ind]
y_coords = y_coords[ind]
z_coords = z_coords[ind]

ind = np.where(y_coords > 0.0)[0]
x_coords = x_coords[ind]
y_coords = y_coords[ind]
z_coords = z_coords[ind]

ind = np.where(z_coords < boxSize)[0]
x_coords = x_coords[ind]
y_coords = y_coords[ind]
z_coords = z_coords[ind]

ind = np.where(z_coords > 0.0)[0]
x_coords = x_coords[ind]
y_coords = y_coords[ind]
z_coords = z_coords[ind]

# count number of particles
N = x_coords.size

# make the coords and radius arrays again
coords_2 = np.zeros((N, 3))
coords_2[:, 0] = x_coords
coords_2[:, 1] = y_coords
coords_2[:, 2] = z_coords

radius = np.sqrt(coords_2[:, 0] ** 2 + coords_2[:, 1] ** 2 + coords_2[:, 2] ** 2)

# Create the file
output_dir = sys.argv[2] 
filename = "Hydrostatic.hdf5"
file_path = os.path.join(output_dir, filename)
file = h5py.File(file_path, "w")

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = boxSize
grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFileOutputsPerSnapshot"] = 1
# grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3

# Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = const_unit_length_in_cgs
grp.attrs["Unit mass in cgs (U_M)"] = const_unit_mass_in_cgs
grp.attrs["Unit time in cgs (U_t)"] = (
    const_unit_length_in_cgs / const_unit_velocity_in_cgs
)
grp.attrs["Unit current in cgs (U_I)"] = 1.0
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0


# Particle group
grp = file.create_group("/PartType0")
grp.create_dataset("Coordinates", data=coords_2, dtype="d")

# All velocities set to zero
v = np.zeros((N, 3))
grp.create_dataset("Velocities", data=v, dtype="f")

# All particles of equal mass
mass = 1.0 / N
m = np.full((N,), mass)
grp.create_dataset("Masses", data=m, dtype="f")

# Smoothing lengths
smoothing_lengths = (4.0 * np.pi * radius ** 2 / N) ** (
    1.0 / 3.0
)  # local mean inter-particle separation
h = np.full((N,), eta * smoothing_lengths)
grp.create_dataset("SmoothingLength", data=h, dtype="f")

# Internal energies
u = v_c ** 2 / (2.0 * (gamma - 1.0))
u = np.full((N,), u)
grp.create_dataset("InternalEnergy", data=u, dtype="f")

# Particle IDs
ids = 1 + np.linspace(0, N, N, endpoint=False)
grp.create_dataset("ParticleIDs", data=ids, dtype="L")

# Material IDs
mat_ids = np.zeros(N)
grp.create_dataset("MaterialIDs", data=mat_ids, dtype="L")

# Density
# NOTE: extract density from 0th snapshot if sys.argv != None
if str(sys.argv[3]) != "None":

    initial_density_folder = sys.argv[3]
    filename = "Hydrostatic_0000.hdf5"
    snapshot_path_initial_density = os.path.join(initial_density_folder, filename)
    with h5py.File(snapshot_path_initial_density, "r") as f:
        rho_init = f["/PartType0/Densities"][:]

    grp.create_dataset("Density", data=rho_init, dtype="f")

else:
    def density_profile(r):
        # Taken from the density profile script
        return r ** -2 / (4.0 * np.pi)
    
    # Determine the densities of all particles based on the radius 
    densities = density_profile(radius)
    grp.create_dataset("Density", data=densities, dtype="f")

file.close()