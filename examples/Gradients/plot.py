################################################################################
# This file is part of SWIFT.
# Copyright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
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

import scipy as sp
import pylab as pl
import numpy as np
import h5py
import sys

# this file plots the gradients of the density in the x and y direction for
# the given input file and saves the result as gradiens_NAME.png

inputfile = sys.argv[1]
outputfile = "gradients_{0}.png".format(sys.argv[2])

f = h5py.File(inputfile, "r")
rho = np.array(f["/PartType0/Density"])
gradrho = np.array(f["/PartType0/GradDensity"])
coords = np.array(f["/PartType0/Coordinates"])

fig, ax = pl.subplots(1,2, sharey=True)

ax[0].plot(coords[:,0], rho, "r.", label="density")
ax[0].plot(coords[:,0], gradrho[:,0], "b.", label="grad density x")
ax[0].set_xlabel("x")
ax[0].legend(loc="best")

ax[1].plot(coords[:,1], rho, "r.", label="density")
ax[1].plot(coords[:,1], gradrho[:,1], "b.", label="grad density y")
ax[1].set_xlabel("y")
ax[1].legend(loc="best")

pl.tight_layout()
pl.savefig(outputfile)
