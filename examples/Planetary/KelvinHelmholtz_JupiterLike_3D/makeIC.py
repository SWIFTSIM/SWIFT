###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2025 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
#               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

# Generates a swift IC file for the Kelvin-Helmholtz test in a periodic box

# Constants
R_earth = 6371000  # Earth radius
R_jupiter = 11.2089 * R_earth  # Jupiter radius

# Parameters
N2_l = 128  # Particles along one edge in the low-density region
N2_depth = 18  # Particles in z direction in low-density region
matID1 = 304  # Central region material ID: AQUA
matID2 = 307  # Outskirts material ID: CD21 H--He
P1 = 3.2e12  # Central region pressure
P2 = 3.2e12  # Outskirts pressure
u1 = 283591514  # Central region specific internal energy
u2 = 804943158  # Outskirts specific internal energy
rho1_approx = 9000  # Central region density. Readjusted later
rho2 = 3500  # Outskirts density
boxsize_l = R_jupiter  # size of simulation box in x and y dimension
v1 = boxsize_l / 10000  # Central region velocity
v2 = -boxsize_l / 10000  # Outskirts velocity
boxsize_depth = boxsize_l * N2_depth / N2_l  # size of simulation box in z dimension
mass = rho2 * (boxsize_l * boxsize_l * boxsize_depth) / (N2_l * N2_l * N2_depth)
fileOutputName = "kelvin_helmholtz.hdf5"
# ---------------------------------------------------

# Start by calculating N1_l and rho1
numPart2 = N2_l * N2_l * N2_depth
numPart1_approx = int(numPart2 / rho2 * rho1_approx)

# Consider numPart1 = N1_l * N1_l * N1_depth
# Substituting boxsize_depth / boxsize_l = N1_depth / N1_l gives,
# numPart1 = N1_l * N1_l * (boxsize_depth / boxsize_l) * N1_l, which ranges to:
N1_l = int(np.cbrt(numPart1_approx * boxsize_l / boxsize_depth))
# Make sure this is a multiple of 4 since this is the number of KH vortices
N1_l -= N1_l % 4

N1_depth = int(boxsize_depth * N1_l / boxsize_l)
numPart1 = int(N1_l * N1_l * N1_depth)

# The density of the central region can then be calculated
rho1 = mass * (N1_l * N1_l * N1_depth) / (boxsize_l * boxsize_l * boxsize_depth)

# Now construct two lattices of particles in the two regions
A2_coords1 = np.empty((numPart1, 3))
A2_coords2 = np.empty((numPart2, 3))
A2_vel1 = np.zeros((numPart1, 3))
A2_vel2 = np.zeros((numPart2, 3))
A2_vel1[:, 0] = v1
A2_vel2[:, 0] = v2

A1_mat1 = np.full(numPart1, matID1)
A1_mat2 = np.full(numPart2, matID2)
A1_m1 = np.full(numPart1, mass)
A1_m2 = np.full(numPart2, mass)
A1_rho1 = np.full(numPart1, rho1)
A1_rho2 = np.full(numPart2, rho2)
A1_u1 = np.full(numPart1, u1)
A1_u2 = np.full(numPart2, u2)
A1_h1 = np.full(numPart1, boxsize_l / N1_l)
A1_h2 = np.full(numPart2, boxsize_l / N2_l)

# Particles in the central region
for i in range(N1_depth):
    for j in range(N1_l):
        for k in range(N1_l):
            index = i * N1_l ** 2 + j * N1_l + k
            A2_coords1[index, 0] = (j / float(N1_l) + 1.0 / (2.0 * N1_l)) * boxsize_l
            A2_coords1[index, 1] = (k / float(N1_l) + 1.0 / (2.0 * N1_l)) * boxsize_l
            A2_coords1[index, 2] = (
                i / float(N1_depth) + 1.0 / (2.0 * N1_depth)
            ) * boxsize_depth

# Particles in the outskirts
for i in range(N2_depth):
    for j in range(N2_l):
        for k in range(N2_l):
            index = i * N2_l ** 2 + j * N2_l + k
            A2_coords2[index, 0] = (j / float(N2_l) + 1.0 / (2.0 * N2_l)) * boxsize_l
            A2_coords2[index, 1] = (k / float(N2_l) + 1.0 / (2.0 * N2_l)) * boxsize_l
            A2_coords2[index, 2] = (
                i / float(N2_depth) + 1.0 / (2.0 * N2_depth)
            ) * boxsize_depth


# Masks for the particles to be selected for the outer and inner regions
mask1 = abs(A2_coords1[:, 1] - 0.5 * boxsize_l) < 0.25 * boxsize_l
mask2 = abs(A2_coords2[:, 1] - 0.5 * boxsize_l) > 0.25 * boxsize_l

# The positions of the particles are now selected
# and the placement of the lattices are adjusted to give appropriate interfaces
A2_coords_inside = A2_coords1[mask1, :]
A2_coords_outside = A2_coords2[mask2, :]

# Calculate the separation of particles across the density discontinuity
pcl_separation_1 = np.cbrt(mass / rho1)
pcl_separation_2 = np.cbrt(mass / rho2)
boundary_separation = 0.5 * (pcl_separation_1 + pcl_separation_2)

# Shift all the "inside" particles to get boundary_separation across the bottom interface
min_y_inside = np.min(A2_coords_inside[:, 1])
max_y_outside_bot = np.max(
    A2_coords_outside[A2_coords_outside[:, 1] - 0.5 * boxsize_l < -0.25 * boxsize_l, 1]
)
shift_distance_bot = boundary_separation - (min_y_inside - max_y_outside_bot)
A2_coords_inside[:, 1] += shift_distance_bot

# Shift the top section of the "outside" particles to get boundary_separation across the top interface
max_y_inside = np.max(A2_coords_inside[:, 1])
min_y_outside_top = np.min(
    A2_coords_outside[A2_coords_outside[:, 1] - 0.5 * boxsize_l > 0.25 * boxsize_l, 1]
)
shift_distance_top = boundary_separation - (min_y_outside_top - max_y_inside)
A2_coords_outside[
    A2_coords_outside[:, 1] - 0.5 * boxsize_l > 0.25 * boxsize_l, 1
] += shift_distance_top

# Adjust box size in y direction based on the shifting of the lattices.
new_box_y = boxsize_l + shift_distance_top

# Now the two lattices can be combined
A2_coords = np.append(A2_coords_inside, A2_coords_outside, axis=0)
A2_vel = np.append(A2_vel1[mask1], A2_vel2[mask2], axis=0)
A1_mat = np.append(A1_mat1[mask1], A1_mat2[mask2], axis=0)
A1_m = np.append(A1_m1[mask1], A1_m2[mask2], axis=0)
A1_rho = np.append(A1_rho1[mask1], A1_rho2[mask2], axis=0)
A1_u = np.append(A1_u1[mask1], A1_u2[mask2], axis=0)
A1_h = np.append(A1_h1[mask1], A1_h2[mask2], axis=0)
numPart = np.size(A1_m)
A1_ids = np.linspace(1, numPart, numPart)

# Finally add the velocity perturbation
vel_perturb_factor = 0.01 * (v1 - v2)
A2_vel[:, 1] = vel_perturb_factor * np.sin(
    2 * np.pi * A2_coords[:, 0] / (0.5 * boxsize_l)
)

# Write ICs to file
with h5py.File(fileOutputName, "w") as f:
    # Header
    grp = f.create_group("/Header")
    grp.attrs["BoxSize"] = [boxsize_l, new_box_y, boxsize_depth]
    grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
    grp.attrs["Time"] = 0.0
    grp.attrs["NumFileOutputsPerSnapshot"] = 1
    grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["Dimension"] = 3

    # Units
    grp = f.create_group("/Units")
    grp.attrs["Unit length in cgs (U_L)"] = 100.0
    grp.attrs["Unit mass in cgs (U_M)"] = 1000.0
    grp.attrs["Unit time in cgs (U_t)"] = 1.0
    grp.attrs["Unit current in cgs (U_I)"] = 1.0
    grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

    # Particle group
    grp = f.create_group("/PartType0")
    ds = grp.create_dataset("Coordinates", (numPart, 3), "d")
    ds[()] = A2_coords
    ds = grp.create_dataset("Velocities", (numPart, 3), "f")
    ds[()] = A2_vel
    ds = grp.create_dataset("Masses", (numPart, 1), "f")
    ds[()] = A1_m.reshape((numPart, 1))
    ds = grp.create_dataset("Density", (numPart, 1), "f")
    ds[()] = A1_rho.reshape((numPart, 1))
    ds = grp.create_dataset("SmoothingLength", (numPart, 1), "f")
    ds[()] = A1_h.reshape((numPart, 1))
    ds = grp.create_dataset("InternalEnergy", (numPart, 1), "f")
    ds[()] = A1_u.reshape((numPart, 1))
    ds = grp.create_dataset("ParticleIDs", (numPart, 1), "L")
    ds[()] = A1_ids.reshape((numPart, 1))
    ds = grp.create_dataset("MaterialIDs", (numPart, 1), "i")
    ds[()] = A1_mat.reshape((numPart, 1))
