################################################################################
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
# along with this program. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

##
# This script plots the Disc-Patch_*.hdf5 snapshots.
# It takes two (optional) parameters: the counter value of the first and last
# snapshot to plot (default: 0 21).
##

import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import pylab as pl
import glob
import sys

# Parameters
surface_density = 10.
scale_height = 100.
z_disc = 400.
z_trunc = 300.
z_max = 350.
utherm = 20.2678457288
gamma = 5. / 3.

start = 0
stop = 21
if len(sys.argv) > 1:
  start = int(sys.argv[1])
if len(sys.argv) > 2:
  stop = int(sys.argv[2])

# Get the analytic solution for the density
def get_analytic_density(x):
  return 0.5 * surface_density / scale_height / \
           np.cosh( (x - z_disc) / scale_height )**2

# Get the analytic solution for the (isothermal) pressure
def get_analytic_pressure(x):
  return (gamma - 1.) * utherm * get_analytic_density(x)

# Get the data fields to plot from the snapshot file with the given name:
#  snapshot time, z-coord, density, pressure, velocity norm
def get_data(name):
  file = h5py.File(name, "r")
  coords = np.array(file["/PartType0/Coordinates"])
  rho = np.array(file["/PartType0/Density"])
  u = np.array(file["/PartType0/InternalEnergy"])
  v = np.array(file["/PartType0/Velocities"])

  P = (gamma - 1.) * rho * u

  vtot = np.sqrt( v[:,0]**2 + v[:,1]**2 + v[:,2]**2 )

  return float(file["/Header"].attrs["Time"]), coords[:,2], rho, P, vtot

# scan the folder for snapshot files and plot all of them (within the requested
# range)
for f in sorted(glob.glob("Disc-Patch_*.hdf5")):
  num = int(f[-8:-5])
  if num < start or num > stop:
    continue

  print "processing", f, "..."

  zrange = np.linspace(0., 2. * z_disc, 1000)
  time, z, rho, P, v = get_data(f)

  fig, ax = pl.subplots(3, 1, sharex = True)

  ax[0].plot(z, rho, "r.")
  ax[0].plot(zrange, get_analytic_density(zrange), "k-")
  ax[0].plot([z_disc - z_max, z_disc - z_max], [0, 10], "k--", alpha=0.5)
  ax[0].plot([z_disc + z_max, z_disc + z_max], [0, 10], "k--", alpha=0.5)
  ax[0].plot([z_disc - z_trunc, z_disc - z_trunc], [0, 10], "k--", alpha=0.5)
  ax[0].plot([z_disc + z_trunc, z_disc + z_trunc], [0, 10], "k--", alpha=0.5)
  ax[0].set_ylim(0., 1.2 * get_analytic_density(z_disc))
  ax[0].set_ylabel("density")

  ax[1].plot(z, v, "r.")
  ax[1].plot(zrange, np.zeros(len(zrange)), "k-")
  ax[1].plot([z_disc - z_max, z_disc - z_max], [0, 10], "k--", alpha=0.5)
  ax[1].plot([z_disc + z_max, z_disc + z_max], [0, 10], "k--", alpha=0.5)
  ax[1].plot([z_disc - z_trunc, z_disc - z_trunc], [0, 10], "k--", alpha=0.5)
  ax[1].plot([z_disc + z_trunc, z_disc + z_trunc], [0, 10], "k--", alpha=0.5)
  ax[1].set_ylim(-0.5, 10.)
  ax[1].set_ylabel("velocity norm")

  ax[2].plot(z, P, "r.")
  ax[2].plot(zrange, get_analytic_pressure(zrange), "k-")
  ax[2].plot([z_disc - z_max, z_disc - z_max], [0, 10], "k--", alpha=0.5)
  ax[2].plot([z_disc + z_max, z_disc + z_max], [0, 10], "k--", alpha=0.5)
  ax[2].plot([z_disc - z_trunc, z_disc - z_trunc], [0, 10], "k--", alpha=0.5)
  ax[2].plot([z_disc + z_trunc, z_disc + z_trunc], [0, 10], "k--", alpha=0.5)
  ax[2].set_xlim(0., 2. * z_disc)
  ax[2].set_ylim(0., 1.2 * get_analytic_pressure(z_disc))
  ax[2].set_xlabel("z")
  ax[2].set_ylabel("pressure")

  pl.suptitle("t = {0:.2f}".format(time))

  pl.savefig("{name}.png".format(name = f[:-5]))
  pl.close()
