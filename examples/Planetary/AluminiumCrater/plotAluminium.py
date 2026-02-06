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
"""
Plots the solution of the square test in a smoothed way using SWIFTsimIO's
smoothed plotting.
"""

import sys

import matplotlib

matplotlib.use("Agg")

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

import h5py

from swiftsimio import load
from swiftsimio.visualisation import project_gas_pixel_grid


import matplotlib
matplotlib.use("Agg")
from pylab import *
from scipy import stats


from mpl_toolkits.axes_grid1 import make_axes_locatable

snap = int(sys.argv[1])

sim = h5py.File("snapshots/impact_%04d.hdf5"%snap, "r")

x = sim["/PartType0/Coordinates"][:,0]
y = sim["/PartType0/Coordinates"][:,1]
z = sim["/PartType0/Coordinates"][:,2]
vx = sim["/PartType0/Velocities"][:,0]
vy = sim["/PartType0/Velocities"][:,1]
vz = sim["/PartType0/Velocities"][:,2]
rho = sim["/PartType0/Densities"][:]
P = sim["/PartType0/Pressures"][:]
u = sim["/PartType0/InternalEnergies"][:]
mat_ids = sim["/PartType0/MaterialIDs"][:]
h = sim["/PartType0/SmoothingLengths"][:]
ids = sim["/PartType0/ParticleIDs"][:]

v = np.sqrt(vx**2 + vy**2 + vz**2)


x = np.array(x)
y = np.array(y)
z = np.array(z)
rho = np.array(rho)
P = np.array(P)
u = np.array(u)
h = np.array(h)
mat_ids = np.array(mat_ids)
v = np.array(v)
ids = np.array(ids)

mid = 0.5 * (np.max(y) + np.min(y))

plot_length=0.5e-3
mask = np.logical_and(y > mid - plot_length, y < mid)

x_midz = x[mask]
y_midz = y[mask]
z_midz = z[mask]
rho_midz = rho[mask]
P_midz = P[mask]
u_midz = u[mask]
h_midz = h[mask]
mat_ids_midz = mat_ids[mask]
v_midz = v[mask]
ids_midz = ids[mask]
# Plot the interesting quantities
fig, ax= plt.subplots(1, 1, figsize=(10, 10))

size = 1#0.05

cmap = plt.get_cmap('rainbow')
norm = mpl.colors.Normalize(vmin=np.min(u_midz),vmax=np.max(u_midz))

ax.scatter(x_midz, z_midz, c=u_midz, norm=norm, cmap=cmap,alpha=1,s=size)
ax.set_aspect('equal', 'box')
#ax.set_xticks([])
#ax.set_yticks([])
ax.set_facecolor((0.8,0.8,0.8))
ax.set_xlim((0.3 * np.max(x), 0.7 * np.max(x)))
ax.set_ylim((0.018,0.023))

sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)

cbar = fig.colorbar(sm, cax=cax, orientation='vertical')
cbar.set_label('Density')

plt.savefig("impact_%04d.png"%snap, dpi=200)
