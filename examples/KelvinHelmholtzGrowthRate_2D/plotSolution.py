###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
import matplotlib
matplotlib.use("Agg")
import pylab as pl
import sys

if len(sys.argv) < 2:
  print "No final snapshot number provided!"
  exit()
lastsnap = int(sys.argv[1])

# Read the simulation data
t = np.zeros(lastsnap + 1)
ey = np.zeros(lastsnap + 1)
for i in range(lastsnap + 1):
  file = h5py.File("kelvinHelmholtzGrowthRate_{0:04d}.hdf5".format(i), 'r')
  t_snap = float(file["/Header"].attrs["Time"])
  vy = file["/PartType0/Velocities"][:,1]
  m = file["/PartType0/Masses"][:]

  ey_snap = 0.5 * m * vy**2

  t[i] = t_snap
  ey[i] = ey_snap.sum()

pl.semilogy(t, ey, "k.")
pl.savefig("kelvinHelmholtzGrowthRate.png")
