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
    description="Generates a swift IC file for the Magnetised Cloud collapse"
)

parser.add_argument(
    "-n",
    "--nparts",
    help="""
         Number of particles to be used in the Evrard collapse.
         """,
    required=False,
    default=100000,
)

args = vars(parser.parse_args())

# Parameters
Rcloud    = 1.0 #4.63e16
Lbox      = 4.0 * Rcloud

gamma = 5.0 / 3.0  # Gas adiabatic index

M     = 1.0 #1.99e33  # total mass of the sphere
#u0 = 0.05 / M  # initial thermal energy
Omega = 2 * pi / 100 #(4.7e5 * 3.154e7)
B0    = 3.0e-5 

cs0       = 1.0  #2e4
inv_rho_c = 1e-4 #1e13

volume_cloud  = (4/3) * pi * Rcloud**3
rho_in  = M / volume_cloud
rho_out = rho_in 

P = rho_in * cs0 * cs0 * sqrt(1.0 + (rho_in * inv_rho_c)**(4/3))

fileName = "magnetised_cloud.hdf5"
numPart_in  = int(args["nparts"])

numPart_out = int(0.01 * Lbox**3 * numPart_in / volume_cloud)

numPart = numPart_in + numPart_out

print(numPart_out)

# Position cloud particles
r         = Rcloud * random.random(numPart_in)
phi       = 2.0 * pi * random.random(numPart_in)
cos_theta = 2.0 * random.random(numPart_in) - 1.0

sin_theta = sqrt(1.0 - cos_theta ** 2)
cos_phi   = cos(phi)
sin_phi   = sin(phi)

pos_in = zeros((numPart_in, 3))
pos_in[:, 0] = r * sin_theta * cos_phi
pos_in[:, 1] = r * sin_theta * sin_phi
pos_in[:, 2] = r * cos_theta

# Shift particles to put the sphere in the centre of the box
pos_in += array([0.5 * Lbox, 0.5 * Lbox, 0.5 * Lbox])

# Position diffuse atmosphere particles
pos_out = Lbox * random.random((numPart_out, 3))

# Aggregate all particles
pos = concatenate((pos_in,pos_out), axis=0)
h   = ones(numPart) * 2.0 * Rcloud / numPart ** (1.0 / 3.0)

# Solid rotation for cloud particles
v = zeros((numPart, 3))
v[:numPart_in,0] = - Omega * r * sin_phi 
v[:numPart_in,1] = Omega * r * cos_phi

# Other attributes
ids = linspace(1, numPart, numPart)
m = ones(numPart) * M / numPart

u = ones(numPart)
u[:numPart_in] *= P / ((gamma-1) * rho_in) 
u[numPart_in:] *= P / ((gamma-1) * rho_out)

B = zeros((numPart, 3))
B[:, 2] = B0

# --------------------------------------------------

# File
file = h5py.File(fileName, "w")

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = [Lbox, Lbox, Lbox]
grp.attrs["NumPart_Total"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 3

# Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.0
grp.attrs["Unit mass in cgs (U_M)"] = 1.0
grp.attrs["Unit time in cgs (U_t)"] = 1.0
grp.attrs["Unit current in cgs (U_I)"] = 1.0
grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

# Particle group
grp = file.create_group("/PartType0")
grp.create_dataset("Coordinates", data=pos, dtype="d")
grp.create_dataset("Velocities", data=v, dtype="f")
grp.create_dataset("MagneticFluxDensity", data=B, dtype="f")
grp.create_dataset("Masses", data=m, dtype="f")
grp.create_dataset("SmoothingLength", data=h, dtype="f")
grp.create_dataset("InternalEnergy", data=u, dtype="f")
grp.create_dataset("ParticleIDs", data=ids, dtype="L")

file.close()
