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

import sys
import h5py
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

snap = int(sys.argv[1])

with h5py.File(f"cylinders_{snap:04d}.hdf5", "r") as f:
    coords = f["/PartType0/Coordinates"][:]
    rho = f["/PartType0/Densities"][:]
x = coords[:, 0]
y = coords[:, 1]

fig, ax = plt.subplots(figsize=(5, 5))

scatter = ax.scatter(x, y, c=rho, cmap="rainbow", s=0.1)
plt.colorbar(scatter, ax=ax, label="Density")
ax.set_aspect("equal")
ax.set_xlim(0, 20)
ax.set_ylim(0, 20)
ax.set_facecolor((0.8, 0.8, 0.8))

plt.savefig(f"cylinders_{snap:04d}.png", dpi=200)
