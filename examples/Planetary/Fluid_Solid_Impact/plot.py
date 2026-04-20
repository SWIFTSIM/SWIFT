###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2019 Josh Borrow (joshua.borrow@durham.ac.uk)
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
import numpy as np
import matplotlib.pyplot as plt

snap = int(sys.argv[1])

with h5py.File(f"impact_{snap:04d}.hdf5", "r") as f:
    x = f["/PartType0/Coordinates"][:, 0]
    y = f["/PartType0/Coordinates"][:, 1]
    z = f["/PartType0/Coordinates"][:, 2]
    rho = f["/PartType0/Densities"][:]
    P = f["/PartType0/Pressures"][:]
    mat_ids = f["/PartType0/MaterialIDs"][:]

# Mid-plane slice
mid = 0.5 * (y.max() + y.min())
plot_depth = 5e-3
mask = (y > mid - plot_depth) & (y < mid)

x_m, z_m, mat_m = x[mask], z[mask], mat_ids[mask]

fig, ax = plt.subplots(figsize=(10, 10))
ax.scatter(x_m[mat_m == 500], z_m[mat_m == 500], c="b", s=10)
ax.scatter(x_m[mat_m == 501], z_m[mat_m == 501], c="r", s=10)
ax.set_aspect("equal", "box")
ax.set_facecolor((0.8, 0.8, 0.8))
ax.set_ylim(0, 1e-2)

plt.savefig(f"impact_{snap:04d}.png", dpi=200)