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

import h5py
import argparse as ap
from numpy import *

parser = ap.ArgumentParser(
    description="Generates a swift IC file for the Evrard collapse"
)

parser.add_argument(
    "-n",
    "--nparts",
    help="""
         Number of particles to be used in the Evrard collapse.
         """,
    required=False,
    default=100000
)

args = vars(parser.parse_args())

# Parameters
gamma = 5. / 3.      # Gas adiabatic index
M = 1.  # total mass of the sphere
R = 1.               # radius of the sphere
u0 = 0.05 / M        # initial thermal energy
fileName = "evrard.hdf5" 
numPart = int(args["nparts"])

r = R * sqrt(random.random(numPart))
phi = 2. * pi * random.random(numPart)
cos_theta = 2. * random.random(numPart) - 1.

sin_theta = sqrt(1. - cos_theta**2)
cos_phi = cos(phi)
sin_phi = sin(phi)

pos = zeros((numPart, 3))
pos[:,0] = r * sin_theta * cos_phi
pos[:,1] = r * sin_theta * sin_phi
pos[:,2] = r * cos_theta

# shift particles to put the sphere in the centre of the box
pos += array([50. * R, 50. * R, 50. * R])

h = ones(numPart) * 2. * R / numPart**(1. / 3.)

# Generate extra arrays
v = zeros((numPart, 3))
ids = linspace(1, numPart, numPart)
m = ones(numPart) * M / numPart
u = ones(numPart) * u0

#--------------------------------------------------

#File
file = h5py.File(fileName, 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = [100. * R, 100. * R, 100. * R]
grp.attrs["NumPart_Total"] =  [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 3

#Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.
grp.attrs["Unit mass in cgs (U_M)"] = 1.
grp.attrs["Unit time in cgs (U_t)"] = 1.
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

#Particle group
grp = file.create_group("/PartType0")
grp.create_dataset('Coordinates', data=pos, dtype='d')
grp.create_dataset('Velocities', data=v, dtype='f')
grp.create_dataset('Masses', data=m, dtype='f')
grp.create_dataset('SmoothingLength', data=h, dtype='f')
grp.create_dataset('InternalEnergy', data=u, dtype='f')
grp.create_dataset('ParticleIDs', data=ids, dtype='L')

file.close()
