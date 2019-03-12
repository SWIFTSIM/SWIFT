#!/usr/bin/env python
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
################################################################################
import numpy as np
import h5py
import matplotlib.pyplot as plt
from scipy.integrate import odeint


lengthrun = 2001
numbpar = 5

radius = np.zeros((numbpar, lengthrun))
time = np.zeros(lengthrun)
for i in range(0, lengthrun):
    Data = h5py.File("hernquist_%04d.hdf5" % i, "r")
    header = Data["Header"]
    time[i] = header.attrs["Time"]
    particles = Data["PartType1"]
    positions = particles["Coordinates"]
    radius[:, i] = positions[:, 0] - 200.0

col = ["b", "r", "c", "y", "k"]

for i in range(0, numbpar):
    plt.plot(time, radius[i, :], col[i])
    plt.axhline(np.max(radius[i, :]), color=col[i], linestyle="--")
    plt.axhline(-np.max(radius[i, :]), color=col[i], linestyle="--")


plt.ylabel("Radial distance (kpc)")
plt.xlabel("Simulation time (internal units)")
plt.savefig("radial_infall.png")
plt.close()
