###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
#                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

import sys
import numpy as np
import unyt
import swiftsimio as sw

# Fixed seed so that the initial conditions are reproducible
np.random.seed(1)

# Generates a swift IC file containing a cartesian distribution of DM particles
# with a density of 1

# Parameters
periodic = 1  # 1 For periodic box
boxSize = 1.0
rho = 1.0
L = int(sys.argv[1])  # Number of particles along one axis
fileName = "uniform_DM_box.hdf5"

# ---------------------------------------------------
numPart = L**3
mass = boxSize**3 * rho / numPart

# --------------------------------------------------

# File
w = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    boxsize=sw.cosmo_array(
        [boxSize, boxSize, boxSize],
        unyt.cm,
        comoving=True,
        scale_factor=1.0,
        scale_exponent=1,
    ),
    dimension=3,
)

# Particle group
v = np.zeros((numPart, 3))

m = np.full(numPart, mass)

ids = np.linspace(0, numPart, numPart, endpoint=False)
coords = np.random.rand(numPart, 3) * boxSize

w.dark_matter.coordinates = coords * unyt.cm * w.a.comoving
w.dark_matter.velocities = v * unyt.cm / unyt.s * w.a.comoving**0
w.dark_matter.masses = m * unyt.g * w.a.comoving**0
w.dark_matter.particle_ids = (ids + 1).astype(np.uint64)

w.write(fileName)
