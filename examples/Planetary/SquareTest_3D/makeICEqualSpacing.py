###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2025 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
#               2025 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
#               2019 Josh Borrow (joshua.borrow@durham.ac.uk)
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
"""Make initial conditions for the 3D square test with equal particle spacing."""

import h5py
import numpy as np
import unyt
import swiftsimio as sw

# Parameters
N_l = 40  # Number of particles on one side
gamma = 5.0 / 3.0  # Gas adiabatic index
rho_in = 4.0  # Density of inner region
rho_out = 1.0  # Density of outer region
P_in = 2.5  # Pressure of inner region
P_out = 2.5  # Pressure of outer region
fileOutputName = "square_equal_spacing.hdf5"

vol = 1.0
numPart = N_l * N_l * N_l

# Set particle masses
A1_pos_l = np.arange(0, 1, 1.0 / N_l)
A3_pos_x, A3_pos_y, A3_pos_z = np.meshgrid(A1_pos_l, A1_pos_l, A1_pos_l)
pos = np.zeros((numPart, 3), dtype=float)
pos[:, 0] = A3_pos_x.flatten()
pos[:, 1] = A3_pos_y.flatten()
pos[:, 2] = A3_pos_z.flatten()

# 3d mask
mask_inside = np.logical_and.reduce(
    [
        A3_pos_x < 0.75,
        A3_pos_x > 0.25,
        A3_pos_y < 0.75,
        A3_pos_y > 0.25,
        A3_pos_z < 0.75,
        A3_pos_z > 0.25,
    ]
)

# Set particle masses
mass_in = rho_in * vol / numPart
mass_out = rho_out * vol / numPart
m = np.ones_like(A3_pos_x) * mass_out
m[mask_inside] = mass_in
m = m.flatten()

# Set approximate particle smoothing lengths
h = np.ones_like(m) / N_l

# Set particle specific internal energies
u_in = P_in / ((gamma - 1) * rho_in)
u_out = P_out / ((gamma - 1) * rho_out)
u = np.ones_like(A3_pos_x) * u_out
u[mask_inside] = u_in
u = u.flatten()

# Set particle densities
rho = np.ones_like(A3_pos_x) * rho_out
rho[mask_inside] = rho_in
rho = rho.flatten()

# Set particle velocities
vel = np.zeros_like(pos)

# Set particle IDs
ids = np.arange(1, numPart + 1)

# Set particle material IDs
mat = np.zeros(numPart)

# Write ICs to file
w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    boxsize=sw.cosmo_array(
        [vol, vol, vol], unyt.cm, comoving=True, scale_factor=1.0, scale_exponent=1
    ),
    dimension=3,
)

w.gas.coordinates = pos * unyt.cm * w.a.comoving
w.gas.velocities = vel * unyt.cm / unyt.s * w.a.comoving**0
w.gas.masses = m * unyt.g * w.a.comoving**0
w.gas.smoothing_lengths = h * unyt.cm * w.a.comoving
w.gas.internal_energy = u * unyt.cm**2 / unyt.s**2 * w.a.comoving**0
w.gas.particle_ids = ids.astype(np.uint64)

w.write(fileOutputName)

# Add the extra fields
with h5py.File(fileOutputName, "r+") as f:
    grp = f["PartType0"]
    grp.create_dataset("Density", data=rho, dtype="f")
    grp.create_dataset("MaterialIDs", data=mat, dtype="i")
