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
import matplotlib

matplotlib.use("Agg")

from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import numpy as np
from scipy import stats

from swiftsimio import load
from swiftsimio.visualisation import project_gas_pixel_grid


import matplotlib
matplotlib.use("Agg")
from pylab import *
from scipy import stats
import h5py


snap = int(sys.argv[1])
sim = h5py.File("snapshots/impact_%04d.hdf5"%snap, "r")


x = sim["/PartType0/Coordinates"][:,0]
y = sim["/PartType0/Coordinates"][:,1]
z = sim["/PartType0/Coordinates"][:,2]
vx = sim["/PartType0/Velocities"][:,0]
vy = sim["/PartType0/Velocities"][:,1]
u = sim["/PartType0/InternalEnergies"][:]
P = sim["/PartType0/Pressures"][:]
rho = sim["/PartType0/Densities"][:]
h = sim["/PartType0/SmoothingLengths"][:]
mat_id = sim["/PartType0/MaterialIDs"][:]
ids = sim["/PartType0/ParticleIDs"][:]
D = sim["/PartType0/Damage"][:]
active_flaws = sim["/PartType0/NumActiveFlaws"][:]
local_strain = sim["/PartType0/LocalStrain"][:]


# Plot the interesting quantities
plt.clf()
fig, ax= plt.subplots(1, 1, figsize=(5, 5))
size = 1

# Choose what to plot
func_plotted = D


norm = mpl.colors.Normalize(vmin=np.min(func_plotted),vmax=np.max(func_plotted) )
sm = plt.cm.ScalarMappable(cmap='rainbow', norm=norm)
sm.set_array([])
divider = make_axes_locatable(ax)
cax = divider.append_axes('right', size='5%', pad=0.05)
cbar = fig.colorbar(sm, cax=cax, orientation='vertical')
cbar.set_label('Damage', rotation=270)


mask = z < 0.100

ax.scatter(x[mask], y[mask], c=func_plotted[mask], norm=norm, cmap='rainbow',alpha=1,s=size,linewidths=0)

ax.set_aspect('equal', 'box')
ax.set_facecolor((0.8,0.8,0.8))
ax.set_xlim((0,20/100))
ax.set_ylim((0,20/100))


savefig("impact_%04d.png"%snap, dpi=200)