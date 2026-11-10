###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 SWIFT developers.
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

"""
Generate plots of the 3D dam break: particles projected on the x-z plane and
coloured by density, plus the particle count to check that the tank is closed.
"""

import h5py
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

snaps = [0, 5, 10, 20, 35, 50]

fig, axes = plt.subplots(2, 3, figsize=(15, 9.5), constrained_layout=True)
for ax, snap in zip(axes.flat, snaps):
    with h5py.File(f"dam_break_{snap:04d}.hdf5", "r") as f:
        boxsize = f["Header"].attrs["BoxSize"]
        time = f["Header"].attrs["Time"][0]
        pos = f["PartType0/Coordinates"][:]
        rho = f["PartType0/Densities"][:]

    sc = ax.scatter(
        pos[:, 0], pos[:, 2], c=rho, s=1, cmap="viridis", vmin=980, vmax=1010
    )
    ax.set_xlim(0, boxsize[0])
    ax.set_ylim(0, boxsize[2])
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("z [m]")
    ax.set_title(f"t = {time:.3f} s, {len(rho)} particles")

fig.colorbar(sc, ax=axes, label=r"Density [kg m$^{-3}$]", shrink=0.6)
fig.savefig("dam_break.png", dpi=200)
