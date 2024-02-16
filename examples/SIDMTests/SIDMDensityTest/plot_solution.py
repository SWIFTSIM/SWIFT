###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2024 Camila Correa (camila.correa@cea.fr)
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
import sys
import numpy as np
import matplotlib.pyplot as plt

sim = h5py.File("hydro_0000.hdf5", "r")
gas_densities = sim["/PartType0/Densities"][:]
gas_h = sim["/PartType0/SmoothingLengths"][:]
sim.close()

sim = h5py.File("sidm_0000.hdf5", "r")
dm_densities = sim["/PartType1/Densities"][:]
dm_h = sim["/PartType1/SIDM_search_radius"][:]
sim.close()

# Let's make a one-to-one comparison

############
# Plot parameters
params = {
    "font.size": 11,
    "font.family": "Times",
    "text.usetex": True,
    "figure.figsize": (5, 2.5),
    "figure.subplot.left": 0.1,
    "figure.subplot.right": 0.95,
    "figure.subplot.bottom": 0.16,
    "figure.subplot.top": 0.93,
    "figure.subplot.wspace": 0.35,
    "figure.subplot.hspace": 0.35,
    "lines.markersize": 1,
    "lines.linewidth": 0.5,
    "figure.max_open_warning": 0,
}
plt.rcParams.update(params)

fig = plt.figure()
ax = plt.subplot(1, 2, 1)
plt.grid("True")

plt.plot(gas_densities, dm_densities, 'o')
x_range = np.arange(np.min(np.log10(gas_densities)),np.max(np.log10(gas_densities))+0.5,0.5)
x_range = 10**x_range
plt.plot(x_range,x_range,'--',color='black')

plt.axis([0, 4e7, 0, 4e7])
plt.xlabel("Gas densities")
plt.ylabel("DM densities")
ax.tick_params(direction='in', axis='both', which='both', pad=4.5)

# #######
ax = plt.subplot(1, 2, 2)
plt.grid("True")

plt.plot(gas_h, dm_h, 'o')
x_range = np.arange(0, 0.21, 0.01)
plt.plot(x_range,x_range,'--',color='black')

plt.axis([0, 0.2, 0, 0.2])
plt.xlabel("Gas $h$")
plt.ylabel("DM $h$")
ax.tick_params(direction='in', axis='both', which='both', pad=4.5)

plt.savefig("./Density_test.png", dpi=300)
plt.close()
############