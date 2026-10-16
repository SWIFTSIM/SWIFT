###############################################################################
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
###############################################################################
"""Uniform lattice of neutral gas at n_H = 1 cm^-3 and 1000 K."""

import h5py
import numpy as np

n_side = 16
box_size = 1.0  # kpc
unit_length = 3.0857e21  # cm
unit_mass = 2.0e33  # g
unit_time = 3.0857e16  # s
m_h = 1.67262171e-24  # g
k_b = 1.3806504e-16  # erg / K
x_h = 0.76
n_h = 1.0  # cm^-3
temperature = 1000.0  # K
gamma = 5.0 / 3.0
mu = 4.0 / (1.0 + 3.0 * x_h)

rho = n_h * m_h / x_h / (unit_mass / unit_length**3)
n = n_side**3
q = (np.arange(n_side) + 0.5) / n_side * box_size
pos = np.array(np.meshgrid(q, q, q, indexing="ij")).reshape(3, -1).T
u = k_b * temperature / (mu * m_h * (gamma - 1.0)) * (unit_time / unit_length) ** 2

with h5py.File("cosmicRayBox.hdf5", "w") as f:
    g = f.create_group("Header")
    g.attrs["BoxSize"] = box_size
    g.attrs["NumPart_Total"] = [n, 0, 0, 0, 0, 0]
    g.attrs["NumPart_Total_HighWord"] = [0] * 6
    g.attrs["NumPart_ThisFile"] = [n, 0, 0, 0, 0, 0]
    g.attrs["Time"] = 0.0
    g.attrs["NumFilesPerSnapshot"] = 1
    g.attrs["MassTable"] = [0.0] * 6
    g.attrs["Flag_Entropy_ICs"] = 0
    g = f.create_group("Units")
    g.attrs["Unit length in cgs (U_L)"] = unit_length
    g.attrs["Unit mass in cgs (U_M)"] = unit_mass
    g.attrs["Unit time in cgs (U_t)"] = unit_time
    g.attrs["Unit current in cgs (U_I)"] = 1.0
    g.attrs["Unit temperature in cgs (U_T)"] = 1.0
    g = f.create_group("PartType0")
    g.create_dataset("Coordinates", data=pos)
    g.create_dataset("Velocities", data=np.zeros((n, 3), dtype="f"))
    g.create_dataset("Masses", data=np.full(n, rho * box_size**3 / n, dtype="f"))
    g.create_dataset("SmoothingLength", data=np.full(n, 1.2348 * box_size / n_side, dtype="f"))
    g.create_dataset("InternalEnergy", data=np.full(n, u, dtype="f"))
    g.create_dataset("ParticleIDs", data=np.arange(1, n + 1, dtype="L"))
