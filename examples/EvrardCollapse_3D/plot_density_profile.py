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
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

# plots the density, velocity and pressure as a function of radius for all
# snapshots in the current directory

import numpy as np
import h5py
import matplotlib
matplotlib.use("Agg")
import pylab as pl
import glob
import os
import sys
import multiprocessing as mp

pl.rcParams["figure.figsize"] = (12, 6)
pl.rcParams["text.usetex"] = True

def plot_file(filename):
  imagefile = "{0}.png".format(filename[:-5])
  if os.path.exists(imagefile) and \
     os.path.getmtime(imagefile) > os.path.getmtime(filename):
    return filename
  file = h5py.File(filename, 'r')
  coords = np.array(file["/PartType0/Coordinates"])
  box = np.array(file["/Header"].attrs["BoxSize"])
  rho = np.array(file["/PartType0/Density"])
  v = np.array(file["/PartType0/Velocities"])
  P = np.array(file["/PartType0/Pressure"])

  coords -= 0.5 * box

  r = np.sqrt(coords[:,0]**2 + coords[:,1]**2 + coords[:,2]**2)
  vr = (coords[:,0] * v[:,0] + coords[:,1] * v[:,1] + coords[:,2] * v[:,2]) / r

  fig, ax = pl.subplots(1, 3, sharex = True)

  ax[0].semilogy(r, rho, "k.")
  ax[1].plot(r, vr, "k.")
  ax[2].semilogy(r, P, "k.")

  for i in range(3):
    ax[i].set_xlabel("radius")
  ax[0].set_ylabel("density")
  ax[1].set_ylabel("velocity")
  ax[2].set_ylabel("pressure")
  ax[1].set_title("time = {0:.2f}".format(
                    	float(file["/Header"].attrs["Time"])))	

  pl.tight_layout()
  pl.savefig(imagefile)
  pl.close()

  return filename

num_threads = 1
first_snap = 0
last_snap = len(glob.glob("evrard_????.hdf5"))

if len(sys.argv) > 1:
  num_threads = int(sys.argv[1])
if len(sys.argv) > 2:
  first_snap = int(sys.argv[2])
if len(sys.argv) > 3:
  last_snap = int(sys.argv[3]) + 1

print "Plotting snapshots {0} to {1} using {2} threads...".format(
  first_snap, last_snap - 1, num_threads)

pool = mp.Pool(num_threads)
results = []
for i in range(first_snap, last_snap):
  filename = "evrard_{0:04d}.hdf5".format(i)
  results.append(pool.apply_async(plot_file, (filename,)))

for result in results:
  print result.get(), "done"
