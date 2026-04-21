###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

# Generates a swift IC file for
# See Antoci+ 2007 (and Gray+ 2001)

L = 0.2  # Length of plate (m)
H = 0.02  # Height of plate (m)
rho_0 = 1000  # Density (kg / m^3)
mat_id = 500  # Material ID
K = 3.25e6  # Bulk modulus (N / m^2)
mu = 7.15e5  # Shear modulus (N / m^2)
c_0 = np.sqrt(K / rho_0)  # Sound speed
kL = 1.875  # Fundamental mode
k = kL / L  # Wavenumber of fundamental mode
v_L0 = 0.01  # Initial velocity of the free end in units of c_0
N_H = 10  # Number of particles in the height direction
N_depth = 15  # Number of particles in the depth direction
assert N_H * L / H == int(N_H * L / H)
N_L = int(N_H * L / H)
N_L_boundary = 5  # Number of layers of boundary particles
boxsize_xy = 2 * L  # box size in x and y directions
boxsize_z = H * (N_depth / N_H)  # box size in z direction
fileOutputName = "plate.hdf5"

# -------------------------------------


def f(x):
    return (np.cos(kL) + np.cosh(kL)) * (np.cosh(k * x) - np.cos(k * x)) + (
        np.sin(kL) - np.sinh(kL)
    ) * (np.sinh(k * x) - np.sin(k * x))


def v_y(x):
    return v_L0 * c_0 * f(x) / f(L)


numPart = (N_L + N_L_boundary) * N_H * N_depth

# Grid
x = np.arange(N_L + N_L_boundary) * (L / N_L)
y = np.arange(N_H) * (H / N_H)
z = np.arange(N_depth) * (boxsize_z / N_depth)
gx, gy, gz = np.meshgrid(x, y, z, indexing="ij")
pos = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()])

# Reposition to centre of box in y
pos[:, 1] -= H / 2
pos[:, 1] += boxsize_xy / 2

# Boundary particle info
N_boundary_particles_print = N_L_boundary * N_H * N_depth
print(f"Boundary particles:{N_boundary_particles_print}")

# Offset to exclude boundary particles for setting velocity
x_0 = (N_L_boundary + 0.5) * L / N_L
vel = np.zeros((numPart, 3))
mask = pos[:, 0] > x_0
vel[mask, 1] = v_y(pos[mask, 0] - x_0)

# Particle quantities
dx = L / N_L
h = np.ones(numPart) * dx * 1.487
m = np.ones(numPart) * rho_0 * dx**3
ids = np.arange(1, numPart + 1)
mat = np.full(numPart, mat_id)
rho = np.full(numPart, rho_0)
u = np.zeros(numPart)

with h5py.File(fileOutputName, "w") as fileOutput:
    grp = fileOutput.create_group("Header")
    grp.attrs["BoxSize"] = [boxsize_xy, boxsize_xy, boxsize_z]
    grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
    grp.attrs["Time"] = 0.0
    grp.attrs["NumFilesPerSnapshot"] = 1
    grp.attrs["MassTable"] = [0.0] * 6
    grp.attrs["Flag_Entropy_ICs"] = [0] * 6
    grp.attrs["Dimension"] = 3

    grp = fileOutput.create_group("Units")
    grp.attrs["Unit length in cgs (U_L)"] = 100.0
    grp.attrs["Unit mass in cgs (U_M)"] = 1000.0
    grp.attrs["Unit time in cgs (U_t)"] = 1.0
    grp.attrs["Unit current in cgs (U_I)"] = 1.0
    grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

    grp = fileOutput.create_group("PartType0")
    grp.create_dataset("Coordinates", data=pos, dtype="d")
    grp.create_dataset("Velocities", data=vel, dtype="f")
    grp.create_dataset("Masses", data=m.reshape(-1, 1), dtype="f")
    grp.create_dataset("Density", data=rho.reshape(-1, 1), dtype="f")
    grp.create_dataset("SmoothingLength", data=h.reshape(-1, 1), dtype="f")
    grp.create_dataset("InternalEnergy", data=u.reshape(-1, 1), dtype="f")
    grp.create_dataset("ParticleIDs", data=ids.reshape(-1, 1), dtype="i")
    grp.create_dataset("MaterialIDs", data=mat.reshape(-1, 1), dtype="i")
