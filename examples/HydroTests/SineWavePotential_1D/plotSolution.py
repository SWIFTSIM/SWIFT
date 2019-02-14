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
##############################################################################

# Plots some quantities for the snapshot file which is passed on as a command
# line argument (full name)

import numpy as np
import h5py
import sys
import matplotlib
matplotlib.use("Agg")
import pylab as pl

pl.rcParams["figure.figsize"] = (12, 10)
pl.rcParams["text.usetex"] = True

# these should be the same as in makeIC.py
uconst = 20.2615290634
cs2 = 2.*uconst/3.
A = 10.

if len(sys.argv) < 2:
  print "Need to provide a filename argument!"
  exit()

fileName = sys.argv[1]

file = h5py.File(fileName, 'r')
coords = np.array(file["/PartType0/Coordinates"])
rho = np.array(file["/PartType0/Density"])
P = np.array(file["/PartType0/Pressure"])
u = np.array(file["/PartType0/InternalEnergy"])
m = np.array(file["/PartType0/Masses"])
vs = np.array(file["/PartType0/Velocities"])
ids = np.array(file["/PartType0/ParticleIDs"])

x = np.linspace(0., 1., 1000)
rho_x = 1000.*np.exp(-0.5*A/np.pi/cs2*np.cos(2.*np.pi*x))

a = A * np.sin(2. * np.pi * x)

apart = A * np.sin(2. * np.pi * coords[:,0])
tkin = -0.5 * np.dot(apart, coords[:,0])
print tkin, 0.5 * np.dot(m, vs[:,0]**2)

ids_reverse = np.argsort(ids)

gradP = np.zeros(P.shape)
for i in range(len(P)):
  iself = int(ids[i])
  corr = 0.
  im1 = iself-1
  if im1 < 0:
    im1 = len(P)-1
    corr = 1.
  ip1 = iself+1
  if ip1 == len(P):
    ip1 = 0
    corr = 1.
  idxp1 = ids_reverse[ip1]
  idxm1 = ids_reverse[im1]
  gradP[i] = (P[idxp1] - P[idxm1])/(coords[idxp1,0] - coords[idxm1,0])

fig, ax = pl.subplots(2, 2, sharex = True)

ax[0][0].plot(coords[:,0], rho, ".", label = "gizmo-mfm")
ax[0][0].plot(x, rho_x, "-", label = "stable solution")
ax[0][0].set_ylabel("$\\rho{}$ (code units)")
ax[0][0].legend(loc = "best")

ax[0][1].plot(x, a, "-", label = "$\\nabla{}\\Phi{}$ external")
ax[0][1].plot(coords[:,0], gradP/rho, ".",
              label = "$\\nabla{}P/\\rho{}$ gizmo-mfm")
ax[0][1].set_ylabel("force (code units)")
ax[0][1].legend(loc = "best")

ax[1][0].axhline(y = uconst, label = "isothermal $u$", color = "k",
                 linestyle = "--")
ax[1][0].plot(coords[:,0], u, ".", label = "gizmo-mfm")
ax[1][0].set_ylabel("$u$ (code units)")
ax[1][0].set_xlabel("$x$ (code units)")
ax[1][0].legend(loc = "best")

#ax[1][1].plot(coords[:,0], m, "y.")
#ax[1][1].set_ylabel("$m_i$ (code units)")
#ax[1][1].set_xlabel("$x$ (code units)")

ax[1][1].axhline(y = 0., label = "target", color = "k",
                 linestyle = "--")
ax[1][1].plot(coords[:,0], vs[:,0], ".", label = "gizmo-mfm")
ax[1][1].set_ylabel("$v_x$ (code units)")
ax[1][1].set_xlabel("$x$ (code units)")
ax[1][1].legend(loc = "best")

pl.tight_layout()
pl.savefig("{fileName}.png".format(fileName = fileName[:-5]))
