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
import woma

# Load EoS tables
woma.load_eos_tables(["CD21_HHe", "AQUA"])

# Generates a swift IC file for the Rayleigh-Taylor instability test

# Constants
R_earth = 6371000  # Earth radius
R_jupiter = 11.2089 * R_earth  # Jupiter radius

# Parameters
N2_l = 64  # Number of particles along one edge in lower region
N2_depth = 18  # Number of particles along in z dimension in lower region
matID1 = 304  # Upper region material ID: AQUA
matID2 = 307  # Lower region material ID: CD21 H--He
rho1_approx = 9000  # Approximate density of upper region. To be recalculated
rho2 = 3500  # Density of lower region
g = -31.44  # Constant gravitational acceleration
P0 = 3.2e12  # Pressure at interface
boxsize_factor = 0.1 * R_jupiter
dv = 0.00025 * boxsize_factor  # Size of velocity perturbation
boxsize_xy = [
    0.5 * boxsize_factor,
    1.0 * boxsize_factor,
]  # Size of the box in x and y dimensions
boxsize_depth = boxsize_xy[0] * N2_depth / N2_l  # Size of simulation box in z dimension
fixed_region = [
    0.05 * boxsize_factor,
    0.95 * boxsize_factor,
]  # y-range of non fixed_region particles
perturbation_region = [
    0.3 * boxsize_factor,
    0.7 * boxsize_factor,
]  # y-range for the velocity perturbation
fileOutputName = "rayleigh_taylor.hdf5"
# ---------------------------------------------------

# Start by generating grids of particles of the two densities
numPart2 = N2_l * N2_l * N2_depth
numPart1 = int(numPart2 / rho2 * rho1_approx)
N1_l = int(np.cbrt(boxsize_xy[0] * numPart1 / boxsize_depth))
N1_l -= N1_l % 4  # Make RT symmetric across centre of both instability regions
N1_depth = int(boxsize_depth * N1_l / boxsize_xy[0])
numPart1 = int(N1_l * N1_l * N1_depth)
numPart = numPart2 + numPart1

# Calculate particle masses and rho1
part_volume_l = (boxsize_xy[0] * 0.5 * boxsize_xy[1] * boxsize_depth) / numPart2
mass = rho2 * part_volume_l
part_volume_h = (boxsize_xy[0] * 0.5 * boxsize_xy[1] * boxsize_depth) / numPart1
rho1 = mass / part_volume_h

# Now construct two lattices of particles in the two regions
A2_coords1 = np.empty((numPart1, 3))
A2_coords2 = np.empty((numPart2, 3))
A2_vel1 = np.zeros((numPart1, 3))
A2_vel2 = np.zeros((numPart2, 3))
A1_mat1 = np.full(numPart1, matID1)
A1_mat2 = np.full(numPart2, matID2)
A1_m1 = np.full(numPart1, mass)
A1_m2 = np.full(numPart2, mass)
A1_rho1 = np.full(numPart1, rho1)
A1_rho2 = np.full(numPart2, rho2)
A1_u1 = np.empty(numPart1)
A1_u2 = np.empty(numPart2)
A1_h1 = np.full(numPart1, boxsize_xy[0] / N1_l)
A1_h2 = np.full(numPart2, boxsize_xy[0] / N2_l)
A1_ids = np.zeros(numPart)

# Set up boundary particle counter
# Boundary particles are set by the N lowest ids of particles, where N is set when configuring swift
boundary_particles = 1

# Particles in the upper region
for i in range(N1_depth):
    for j in range(N1_l):
        for k in range(N1_l):
            index = i * N1_l ** 2 + j * N1_l + k
            A2_coords1[index, 0] = (j / float(N1_l) + 1.0 / (2.0 * N1_l)) * boxsize_xy[
                0
            ]
            A2_coords1[index, 1] = (k / float(N1_l) + 1.0 / (2.0 * N1_l)) * (
                0.5 * boxsize_xy[1]
            ) + 0.5 * boxsize_xy[1]
            A2_coords1[index, 2] = (
                i / float(N1_depth) + 1.0 / (2.0 * N1_depth)
            ) * boxsize_depth
            A1_rho1[index] = rho1

            # If in top and bottom where particles are fixed
            if (
                A2_coords1[index, 1] < fixed_region[0]
                or A2_coords1[index, 1] > fixed_region[1]
            ):
                A1_ids[index] = boundary_particles
                boundary_particles += 1


# Particles in the lower region
for i in range(N2_depth):
    for j in range(N2_l):
        for k in range(N2_l):
            index = i * N2_l ** 2 + j * N2_l + k
            A2_coords2[index, 0] = (j / float(N2_l) + 1.0 / (2.0 * N2_l)) * boxsize_xy[
                0
            ]
            A2_coords2[index, 1] = (k / float(N2_l) + 1.0 / (2.0 * N2_l)) * (
                0.5 * boxsize_xy[1]
            )
            A2_coords2[index, 2] = (
                i / float(N2_depth) + 1.0 / (2.0 * N2_depth)
            ) * boxsize_depth
            A1_rho2[index] = rho2

            # If in top and bottom where particles are fixed
            if (
                A2_coords2[index, 1] < fixed_region[0]
                or A2_coords2[index, 1] > fixed_region[1]
            ):
                A1_ids[index + numPart1] = boundary_particles
                boundary_particles += 1

print(
    "You need to compile the code with "
    "--enable-boundary-particles=%i" % boundary_particles
)

# Set IDs of non-boundary particles
A1_ids[A1_ids == 0] = np.linspace(
    boundary_particles, numPart, numPart - boundary_particles + 1
)

# The placement of the lattices are now adjusted to give appropriate interfaces
# Calculate the separation of particles across the density discontinuity
pcl_separation_2 = np.cbrt(mass / rho2)
pcl_separation_1 = np.cbrt(mass / rho1)
boundary_separation = 0.5 * (pcl_separation_2 + pcl_separation_1)

# Shift top lattice
min_y1 = min(A2_coords1[:, 1])
max_y2 = max(A2_coords2[:, 1])
shift_distance = boundary_separation - (min_y1 - max_y2)
A2_coords1[:, 1] += shift_distance

# Calculate internal energies
A1_P1 = P0 + g * A1_rho1 * (A2_coords1[:, 1] - 0.5 * boxsize_xy[1])
A1_P2 = P0 + g * A1_rho2 * (A2_coords2[:, 1] - 0.5 * boxsize_xy[1])
A1_u1 = woma.A1_Z_rho_Y(A1_rho1, A1_P1, A1_mat1, Z_choice="u", Y_choice="P")
A1_u2 = woma.A1_Z_rho_Y(A1_rho2, A1_P2, A1_mat2, Z_choice="u", Y_choice="P")

# Now the two lattices can be combined
A2_coords = np.append(A2_coords1, A2_coords2, axis=0)
A2_vel = np.append(A2_vel1, A2_vel2, axis=0)
A1_mat = np.append(A1_mat1, A1_mat2, axis=0)
A1_m = np.append(A1_m1, A1_m2, axis=0)
A1_rho = np.append(A1_rho1, A1_rho2, axis=0)
A1_u = np.append(A1_u1, A1_u2, axis=0)
A1_h = np.append(A1_h1, A1_h2, axis=0)

# Add velocity perturbation
mask_perturb = np.logical_and(
    A2_coords[:, 1] > perturbation_region[0], A2_coords[:, 1] < perturbation_region[1]
)
A2_vel[mask_perturb, 1] = (
    dv
    * (1 + np.cos(8 * np.pi * (A2_coords[mask_perturb, 0] / (boxsize_factor) + 0.25)))
    * (1 + np.cos(5 * np.pi * (A2_coords[mask_perturb, 1] / (boxsize_factor) - 0.5)))
)


# Write ICs to file
with h5py.File(fileOutputName, "w") as f:
    # Header
    grp = f.create_group("/Header")
    grp.attrs["BoxSize"] = [
        boxsize_xy[0],
        boxsize_xy[1] + shift_distance,
        boxsize_depth,
    ]
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
