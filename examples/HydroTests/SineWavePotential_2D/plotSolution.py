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
import pylab as pl

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
u = np.array(file["/PartType0/InternalEnergy"])
agrav = np.array(file["/PartType0/GravAcceleration"])
m = np.array(file["/PartType0/Masses"])
ids = np.array(file["/PartType0/ParticleIDs"])

# ids_reverse gives the index original particle 0 now has in the particle arrays
# and so on
# note that this will only work if the particles do not move away too much from
# there initial positions
ids_reverse = np.argsort(ids)

x = np.linspace(0., 1., 1000)
rho_x = 1000.*np.exp(-0.5*A/np.pi/cs2*np.cos(2.*np.pi*x))

P = cs2*rho

n1D = int(np.sqrt(len(P)))
gradP = np.zeros(P.shape)
for i in range(len(P)):
  iself = int(ids[i]/n1D)
  jself = int(ids[i]-n1D*iself)
  corr = 0.
  im1 = iself-1
  if im1 < 0:
    im1 = n1D-1
    corr = 1.
  ip1 = iself+1
  if ip1 == n1D:
    ip1 = 0
    corr = 1.
  idxp1 = ids_reverse[ip1*n1D+jself]
  idxm1 = ids_reverse[im1*n1D+jself]
  gradP[i] = (P[idxp1]-P[idxm1])/(coords[idxp1,0]-coords[idxm1,0]+corr)

fig, ax = pl.subplots(2, 2)

ax[0][0].plot(coords[:,0], rho, "r.", markersize = 0.5)
ax[0][0].plot(x, rho_x, "g-")
ax[0][1].plot(coords[:,0], gradP/rho, "b.")
ax[1][0].plot(coords[:,0], agrav[:,0], "g.", markersize = 0.5)
ax[1][1].plot(coords[:,0], m, "y.")
pl.savefig("{fileName}.png".format(fileName = fileName[:-5]))
