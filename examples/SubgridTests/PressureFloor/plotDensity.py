#!/usr/bin/env python
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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

import matplotlib.pyplot as plt
from h5py import File
import sys

N = 100
filename = "pressureFloor_%04i.hdf5" % int(sys.argv[-1])

f = File(filename, "r")

rho = f["PartType0"]["Densities"][:]

plt.hist(rho, 100, log=True)
plt.xlabel("Density [M$_\odot$ kpc$^{-3}$]")
plt.ylabel("Number of particles")
plt.savefig("density.png", dpi=200)
