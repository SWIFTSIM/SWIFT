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
"""Make initial conditions for the 3D square test with equal particle mass."""

import h5py
import numpy as np

# Parameters
N_l = 40  # Number of particles on one side
gamma = 5.0 / 3.0  # Gas adiabatic index
rho_in_approx = 16 ** 3 / 10 ** 3  # Density of inner region
rho_out = 1.0  # Density of outer region
P_in = 2.5  # Pressure of inner region
P_out = 2.5  # Pressure of outer region
fileOutputName = "square_equal_mass.hdf5"

vol = 1.0
cube_vol_factor = 1 / 8
numPart_out = N_l * N_l * N_l
numPart_in_approx = N_l * N_l * N_l * cube_vol_factor * rho_in_approx / rho_out

# Calculate the number of particles on one side of relevant region
N_l_out = int(np.cbrt(numPart_out))
N_l_in = int(np.cbrt(numPart_in_approx))
numPart_in = int(N_l_in ** 3)

# Set up outer region (cube not yet removed from this lattice)
pos_out = np.zeros((numPart_out, 3))
for i in range(N_l_out):
    for j in range(N_l_out):
        for k in range(N_l_out):
            index = i * N_l_out * N_l_out + j * N_l_out + k
            pos_out[index, 0] = i / (float(N_l_out)) + 1.0 / (2.0 * N_l_out)
            pos_out[index, 1] = j / (float(N_l_out)) + 1.0 / (2.0 * N_l_out)
            pos_out[index, 2] = k / (float(N_l_out)) + 1.0 / (2.0 * N_l_out)

h_out = np.ones(numPart_out) * (1.0 / N_l_out)
m_out = np.ones(numPart_out) * vol * rho_out / numPart_out
u_out = np.ones(numPart_out) * P_out / (rho_out * (gamma - 1.0))
rho_out = np.ones(numPart_out) * rho_out

# Set up inner region
rho_in = m_out[0] * numPart_in / cube_vol_factor
pos_in = np.zeros((numPart_in, 3))
for i in range(N_l_in):
    for j in range(N_l_in):
        for k in range(N_l_in):
            index = i * N_l_in * N_l_in + j * N_l_in + k
            pos_in[index, 0] = 0.25 + i / float(2 * N_l_in) + 1.0 / (2.0 * 2 * N_l_in)
            pos_in[index, 1] = 0.25 + j / float(2 * N_l_in) + 1.0 / (2.0 * 2 * N_l_in)
            pos_in[index, 2] = 0.25 + k / float(2 * N_l_in) + 1.0 / (2.0 * 2 * N_l_in)

h_in = np.ones(numPart_in) * (1.0 / N_l_in)
m_in = np.ones(numPart_in) * m_out[0]
u_in = np.ones(numPart_in) * P_in / (rho_in * (gamma - 1.0))
rho_in = np.ones(numPart_in) * rho_in

# Remove the particles within the central cube from the outer region
mask_out = np.logical_or.reduce(
    (
        pos_out[:, 0] < 0.25,
        pos_out[:, 0] > 0.75,
        pos_out[:, 1] < 0.25,
        pos_out[:, 1] > 0.75,
        pos_out[:, 2] < 0.25,
        pos_out[:, 2] > 0.75,
    )
)
pos_out = pos_out[mask_out, :]
h_out = h_out[mask_out]
u_out = u_out[mask_out]
m_out = m_out[mask_out]
rho_out = rho_out[mask_out]

# Combine inner and outer regions
pos = np.append(pos_out, pos_in, axis=0)
h = np.append(h_out, h_in, axis=0)
u = np.append(u_out, u_in)
m = np.append(m_out, m_in)
rho = np.append(rho_out, rho_in)
numPart = np.size(h)
vel = np.zeros((numPart, 3))
ids = np.linspace(1, numPart, numPart)
mat = np.zeros(numPart)

# Write ICs to file
with h5py.File(fileOutputName, "w") as f:
    # Header
    grp = f.create_group("/Header")
    grp.attrs["BoxSize"] = [vol, vol, vol]
    grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
    grp.attrs["Time"] = 0.0
    grp.attrs["NumFilesPerSnapshot"] = 1
    grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["Dimension"] = 3

    # Units
    grp = f.create_group("/Units")
    grp.attrs["Unit length in cgs (U_L)"] = 1.0
    grp.attrs["Unit mass in cgs (U_M)"] = 1.0
    grp.attrs["Unit time in cgs (U_t)"] = 1.0
    grp.attrs["Unit current in cgs (U_I)"] = 1.0
    grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

    # Particle group
    grp = f.create_group("/PartType0")
    grp.create_dataset("Coordinates", data=pos, dtype="d")
    grp.create_dataset("Velocities", data=vel, dtype="f")
    grp.create_dataset("Masses", data=m.reshape((numPart, 1)), dtype="f")
    grp.create_dataset("Density", data=rho.reshape((numPart, 1)), dtype="f")
    grp.create_dataset("SmoothingLength", data=h.reshape((numPart, 1)), dtype="f")
    grp.create_dataset("InternalEnergy", data=u.reshape((numPart, 1)), dtype="f")
    grp.create_dataset("ParticleIDs", data=ids.reshape((numPart, 1)), dtype="L")
    grp.create_dataset("MaterialIDs", data=mat.reshape((numPart, 1)), dtype="i")
