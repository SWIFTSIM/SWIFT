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
import numpy as np
import h5py
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

snap = int(sys.argv[1])

with h5py.File(f"impact_{snap:04d}.hdf5", "r") as sim:
    x = sim["/PartType0/Coordinates"][:, 0]
    y = sim["/PartType0/Coordinates"][:, 1]
    z = sim["/PartType0/Coordinates"][:, 2]
    vx = sim["/PartType0/Velocities"][:, 0]
    vy = sim["/PartType0/Velocities"][:, 1]
    vz = sim["/PartType0/Velocities"][:, 2]
    rho = sim["/PartType0/Densities"][:]
    P = sim["/PartType0/Pressures"][:]
    u = sim["/PartType0/InternalEnergies"][:]
    mat_ids = sim["/PartType0/MaterialIDs"][:]
    h = sim["/PartType0/SmoothingLengths"][:]
    ids = sim["/PartType0/ParticleIDs"][:]

v = np.sqrt(vx**2 + vy**2 + vz**2)

# Mid-plane slice
mid = 0.5 * (y.max() + y.min())
plot_depth = 0.5e-3
mask = (y > mid - plot_depth) & (y < mid)
x_m, y_m, z_m, rho_m = x[mask], y[mask], z[mask], rho[mask]

# Sort by y
order = np.argsort(y_m)
x_m = x_m[order]
y_m = y_m[order]
z_m = z_m[order]
rho_m = rho_m[order]

# Plot
fig, ax = plt.subplots(figsize=(10, 10))

scatter = ax.scatter(
    x_m, z_m, c=rho_m, cmap="rainbow", s=1, vmin=rho_m.min(), vmax=rho_m.max()
)

ax.set_aspect("equal", "box")
ax.set_facecolor((0.8, 0.8, 0.8))
ax.set_xlim(0.2 * x.max(), 0.8 * x.max())
ax.set_ylim(0.007, 0.015)

# Colorbar
divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.05)
cbar = fig.colorbar(scatter, cax=cax)
cbar.set_label("Density")

plt.savefig(f"impact_{snap:04d}.png", dpi=200)
