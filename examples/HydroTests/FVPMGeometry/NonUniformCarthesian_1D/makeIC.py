################################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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
import numpy as np
import argparse


def parse_options():
    """Parse command line arguments for the ICs generator."""
    parser = argparse.ArgumentParser(
        description="Create a dual-resolution advection test for SWIFT."
    )
    parser.add_argument("-l", "--level", type=int, default=5, help="Resolution level.")
    parser.add_argument(
        "-d",
        "--dimension",
        type=int,
        default=1,
        choices=[1, 2, 3],
        help="Dimensionality.",
    )
    parser.add_argument(
        "-o", "--output", type=str, default="advection.hdf5", help="Output filename."
    )
    parser.add_argument(
        "-c",
        "--cartesian",
        action="store_true",
        help="Use a uniform Cartesian grid (no resolution jump).",
    )

    return parser.parse_args()


########################################
# Main
########################################

args = parse_options()

# Define standard units
UnitMass_in_cgs = 1
UnitLength_in_cgs = 1
UnitVelocity_in_cgs = 1
UnitCurrent_in_cgs = 1  # Amperes
UnitTemp_in_cgs = 1  # Kelvin
UnitTime_in_cgs = UnitLength_in_cgs / UnitVelocity_in_cgs

# Grid parameters
level = args.level
dim = args.dimension
nx = 2**level
dx = 1.0 / nx

# Set the physical parameters
rho_high = 2.0  # Density of the square wave
rho_low = 1.0  # Density of the background
velocity_x = 1.0  # Constant advection

###################
# We want a box that is L=1 in Y and Z, but we split X into two zones
# Left: spacing dx, Right: spacing 2*dx
# To keep it simple, let's make the box X-length = 1.0
lx, ly, lz = 1.0, (1.0 if dim > 1 else dx), (1.0 if dim > 2 else dx)

if args.cartesian:
    # Pure uniform Cartesian grid
    x_coords = np.arange(0, 1.0, dx)
else:
    # Dual-resolution zones
    x_left = np.arange(0, 0.5, dx)
    x_right = np.arange(0.5, 1.0, 2 * dx)
    x_coords = np.concatenate([x_left, x_right])

y_coords = np.arange(0, ly, dx) if dim >= 2 else np.array([0.5 * ly])
z_coords = np.arange(0, lz, dx) if dim == 3 else np.array([0.5 * lz])

# Create the grid
X, Y, Z = np.meshgrid(x_coords, y_coords, z_coords, indexing="ij")
pos = np.vstack([X.ravel(), Y.ravel(), Z.ravel()]).T
N = len(pos)
###################

# Identify the "Blob" - let's put it between x = 0.1 and x = 0.4
blob_mask = (pos[:, 0] > 0.1) & (pos[:, 0] < 0.4)

# Mass calculation: m = rho(x) * dV(x)
mass = np.zeros(N)
rho_field = np.ones(N) * rho_low
rho_field[blob_mask] = rho_high

if args.cartesian:
    # Uniform volume dV = dx^dim
    mass = rho_field * (dx**dim)
    h = np.ones(N) * 1.2348 * dx
else:
    # Dual-resolution volume dV
    left_mask = pos[:, 0] < 0.5
    dv = np.zeros(N)
    dv[left_mask] = dx**dim
    dv[~left_mask] = 2.0 * dx * (dx ** (dim - 1))

    mass = rho_field * dv

    # Smoothing length must follow the grid spacing (dV^(1/dim))
    h = np.zeros(N)
    h[left_mask] = 1.2348 * dx
    h[~left_mask] = 1.2348 * (2.0 * dx)

rho = rho_field

# Create the other arrays
vel = np.zeros((N, 3))
vel[:, 0] = velocity_x
ids = np.arange(N, dtype="L")
u = np.ones(N) * 1.0  # Internal energy
rho = np.ones(N) * rho

########################################
# Write HDF5
with h5py.File(args.output, "w") as f:
    # Header
    grp = f.create_group("/Header")
    grp.attrs["BoxSize"] = [lx, ly, lz]
    grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
    grp.attrs["Time"] = 0.0
    grp.attrs["NumFileOutputsPerSnapshot"] = 1
    grp.attrs["MassTable"] = [0.0] * 6
    grp.attrs["Flag_Entropy_ICs"] = [0] * 6
    grp.attrs["Dimension"] = dim

    grp = f.create_group("/Units")
    grp.attrs["Unit length in cgs (U_L)"] = UnitLength_in_cgs
    grp.attrs["Unit mass in cgs (U_M)"] = UnitMass_in_cgs
    grp.attrs["Unit time in cgs (U_t)"] = UnitTime_in_cgs
    grp.attrs["Unit current in cgs (U_I)"] = UnitCurrent_in_cgs
    grp.attrs["Unit temperature in cgs (U_T)"] = UnitTemp_in_cgs

    # Particles
    p_grp = f.create_group("/PartType0")
    p_grp.create_dataset("Coordinates", data=pos, dtype="d")
    p_grp.create_dataset("Velocities", data=vel, dtype="f")
    p_grp.create_dataset("Masses", data=mass, dtype="f")
    p_grp.create_dataset("SmoothingLength", data=h, dtype="f")
    p_grp.create_dataset("InternalEnergy", data=u, dtype="f")
    p_grp.create_dataset("ParticleIDs", data=ids, dtype="L")
    p_grp.create_dataset("Densities", data=rho, dtype="f")

print(f"File {args.output} generated with {N} particles.")
