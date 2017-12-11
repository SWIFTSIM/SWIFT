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
import scipy.stats as stats

pl.rcParams["figure.figsize"] = (12, 10)
pl.rcParams["text.usetex"] = True

force_replot = False
if len(sys.argv) > 2:
  force_replot = sys.argv[2] in ["True", "true", "yes", "sure", "OK"]
gamma = 5. / 3.

global_rholim = [2., 1.e11]
global_vlim = [-2.6, 7.]
global_Plim = [6.e-4, 4.1e10]
global_Alim = [4.e-15, 3.e-4]

def bin_quantity(x, y, xbin_edges):
  ybin, _, _ = stats.binned_statistic(x, y, bins = xbin_edges,
                                      statistic = "mean")
  y2bin, _, _ = stats.binned_statistic(x, y**2, bins = xbin_edges,
                                       statistic = "mean")
  ysigma = np.sqrt(y2bin - ybin**2)
  return ybin, ysigma

def plot_file(filename):
  imagefile = "{0}.png".format(filename[:-5])
  if (not force_replot) and os.path.exists(imagefile) and \
     os.path.getmtime(imagefile) > os.path.getmtime(filename):
    return filename, global_rholim, global_vlim, global_Plim, global_Alim
  file = h5py.File(filename, 'r')
  coords = np.array(file["/PartType0/Coordinates"])
  box = np.array(file["/Header"].attrs["BoxSize"])
  rho = np.array(file["/PartType0/Density"])
  v = np.array(file["/PartType0/Velocities"])
  P = np.array(file["/PartType0/Pressure"])

  coords -= 0.5 * box

  r = np.sqrt(coords[:,0]**2 + coords[:,1]**2 + coords[:,2]**2)
  vr = (coords[:,0] * v[:,0] + coords[:,1] * v[:,1] + coords[:,2] * v[:,2]) / r

  rbin_edges = np.linspace(r.min(), r.max(), 101)
  rbin_mid = 0.5 * (rbin_edges[1:] + rbin_edges[:-1])

  rhobin, rhosigma = bin_quantity(r, rho, rbin_edges)
  vbin, vsigma = bin_quantity(r, vr, rbin_edges)
  Pbin, Psigma = bin_quantity(r, P, rbin_edges)
  Abin, Asigma = bin_quantity(r, P / rho**gamma, rbin_edges)

  fig, ax = pl.subplots(2, 2, sharex = True)

  ax[0][0].loglog(rbin_mid, rhobin, "-")
  ax[0][0].fill_between(rbin_mid, rhobin - rhosigma, rhobin + rhosigma,
                        alpha = 0.5)
  ax[0][1].semilogx(rbin_mid, vbin, "-")
  ax[0][1].fill_between(rbin_mid, vbin - vsigma, vbin + vsigma,
                        alpha = 0.5)
  ax[1][0].loglog(rbin_mid, Pbin, "-")
  ax[1][0].fill_between(rbin_mid, Pbin - Psigma, Pbin + Psigma,
                        alpha = 0.5)
  ax[1][1].loglog(rbin_mid, Abin, "-")
  ax[1][1].fill_between(rbin_mid, Abin - Asigma, Abin + Asigma,
                        alpha = 0.5)

  ax[0][0].set_yscale("log", nonposy = "clip")
  ax[1][0].set_yscale("log", nonposy = "clip")
  ax[1][1].set_yscale("log", nonposy = "clip")

  ax[0][0].set_xlim(1.e-3, 0.5 * box[0])
  ax[0][0].set_ylim(global_rholim[0], global_rholim[1])
  ax[0][1].set_ylim(global_vlim[0], global_vlim[1])
  ax[1][0].set_ylim(global_Plim[0], global_Plim[1])
  ax[1][1].set_ylim(global_Alim[0], global_Alim[1])

  for i in range(2):
    ax[1][i].set_xlabel("radius")
  ax[0][0].set_ylabel("density")
  ax[0][1].set_ylabel("velocity")
  ax[1][0].set_ylabel("pressure")
  ax[1][1].set_ylabel("entropy")
  ax[0][0].set_title("time = {0:.2f}".format(
                       float(file["/Header"].attrs["Time"])))

  pl.tight_layout()
  pl.savefig(imagefile)
  pl.close()

  rhobin = rhobin[np.isfinite(rhosigma)]
  vbin = vbin[np.isfinite(vsigma)]
  Pbin = Pbin[np.isfinite(Psigma)]
  Abin = Abin[np.isfinite(Asigma)]
  rhosigma = rhosigma[np.isfinite(rhosigma)]
  vsigma = vsigma[np.isfinite(vsigma)]
  Psigma = Psigma[np.isfinite(Psigma)]
  Asigma = Asigma[np.isfinite(Asigma)]

  rhomax = (rhobin + rhosigma).max()
  rhomin = np.where(rhobin - rhosigma > 0.,
                    rhobin - rhosigma,
                    np.ones(len(rhobin)) * rhomax).min()
  vmax = (vbin + vsigma).max()
  vmin = (vbin - vsigma).min()
  Pmax = (Pbin + Psigma).max()
  Pmin = np.where(Pbin - Psigma > 0.,
                  Pbin - Psigma,
                  np.ones(len(Pbin)) * Pmax).min()
  Amax = (Abin + Asigma).max()
  Amin = np.where(Abin - Asigma > 0.,
                  Abin - Asigma,
                  np.ones(len(Abin)) * Amax).min()

  return filename, [rhomin, rhomax], [vmin, vmax], [Pmin, Pmax], [Amin, Amax]

num_threads = 1
first_snap = 0
last_snap = len(glob.glob("evrard_????.hdf5"))

if len(sys.argv) > 1:
  num_threads = int(sys.argv[1])
if len(sys.argv) > 3:
  first_snap = int(sys.argv[3])
if len(sys.argv) > 4:
  last_snap = int(sys.argv[4]) + 1

print "Plotting snapshots {0} to {1} using {2} threads...".format(
  first_snap, last_snap - 1, num_threads)

pool = mp.Pool(num_threads)
results = []
for i in range(first_snap, last_snap):
  filename = "evrard_{0:04d}.hdf5".format(i)
  results.append(pool.apply_async(plot_file, (filename,)))

rholim = [np.inf, -np.inf]
vlim = [np.inf, -np.inf]
Plim = [np.inf, -np.inf]
Alim = [np.inf, -np.inf]
for result in results:
  name, rhol, vl, Pl, Al = result.get()
  print name, "done."
  rholim[0] = min(rholim[0], rhol[0])
  rholim[1] = max(rholim[1], rhol[1])
  vlim[0] = min(vlim[0], vl[0])
  vlim[1] = max(vlim[1], vl[1])
  Plim[0] = min(Plim[0], Pl[0])
  Plim[1] = max(Plim[1], Pl[1])
  Alim[0] = min(Alim[0], Al[0])
  Alim[1] = max(Alim[1], Al[1])

print "rholim:", rholim
print "vlim:", vlim
print "Plim:", Plim
print "Alim:", Alim

energy = np.loadtxt("energy.txt")

pl.rcParams["figure.figsize"] = (8, 6)

pl.plot(energy[:, 0], energy[:, 2], label = "Total energy")
pl.plot(energy[:, 0], energy[:, 3], label = "Kinetic energy")
pl.plot(energy[:, 0], energy[:, 4], label = "Internal energy")
pl.plot(energy[:, 0], energy[:, 5], label = "Potential energy")

pl.xlabel("time")
pl.ylabel("energy")
pl.legend(loc = "best")
pl.tight_layout()
pl.savefig("evrard_energy.png")
