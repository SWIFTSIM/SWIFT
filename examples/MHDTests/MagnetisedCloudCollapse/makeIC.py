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
from scipy.spatial.transform import Rotation

parser = ap.ArgumentParser(
    description="Generates a swift IC file for the Magnetised Cloud collapse"
)

parser.add_argument(
    "-n",
    "--nparts",
    help="""
         Number of particles on a side of the cube out of which we cut the dense collapsing sphere.
         """,
    required=False,
    default=64,
)

args = vars(parser.parse_args())

# Parameters
Rcloud = 4.628516371e16
Lbox = 10.0 * Rcloud

gamma = 5.0 / 3.0  # Gas adiabatic index

M = 1.99e33  # total mass of the sphere
Omega = 2 * pi / (4.7e5 * 3.1536e7)

mu_fid = 10.0
Bini_fid = 6.1e-7

mu = 5
Bini = Bini_fid * (mu_fid / mu)

cs0 = 2e4
inv_rho_c = 1e14

volume_cloud = (4 / 3) * pi * Rcloud ** 3
volume_cloud_box = (2 * Rcloud) ** 3
volume_sim_box = Lbox ** 3

rho_in = M / volume_cloud
rho_out = rho_in / 360

P_in = rho_in * cs0 * cs0 * sqrt(1.0 + (rho_in * inv_rho_c) ** (4 / 3))
P_out = rho_out * cs0 * cs0 * sqrt(1.0 + (rho_out * inv_rho_c) ** (4 / 3))

fileName = "magnetised_cloud.hdf5"

numPart_in_side = int(args["nparts"])
numPart_out_side = int(
    floor(numPart_in_side * cbrt(volume_sim_box / volume_cloud_box / 360))
)

# Position cloud particles

x_ = linspace(-Rcloud, Rcloud, num=numPart_in_side, endpoint=True)
y_ = linspace(-Rcloud, Rcloud, num=numPart_in_side, endpoint=True)
z_ = linspace(-Rcloud, Rcloud, num=numPart_in_side, endpoint=True)

x, y, z = meshgrid(x_, y_, z_, indexing="ij")

x = x.flatten()
y = y.flatten()
z = z.flatten()

r = stack((x, y, z), axis=1)
pos_in = r[r[:, 0] ** 2 + r[:, 1] ** 2 + r[:, 2] ** 2 < Rcloud ** 2]

numPart_in = int(pos_in.shape[0])

phi = arctan2(pos_in[:, 1], pos_in[:, 0])
cos_phi = cos(phi)
sin_phi = sin(phi)

R = sqrt(pos_in[:, 0] ** 2 + pos_in[:, 1] ** 2)

# Shift particles to put the sphere in the centre of the box
pos_in += array([0.5 * Lbox, 0.5 * Lbox, 0.5 * Lbox])

# Position diffuse atmosphere particles

step = Lbox / (numPart_out_side - 1)
x_ = arange(0.0, Lbox, step) + step / 2  # , endpoint=True)
y_ = arange(0.0, Lbox, step) + step / 2  # , endpoint=True)
z_ = arange(0.0, Lbox, step) + step / 2  # , endpoint=True)

x, y, z = meshgrid(x_, y_, z_, indexing="ij")

x = x.flatten()
y = y.flatten()
z = z.flatten()

pos_out = stack((x, y, z), axis=1)
pos_out = pos_out[
    (pos_out[:, 0] - 0.5 * Lbox) ** 2
    + (pos_out[:, 1] - 0.5 * Lbox) ** 2
    + (pos_out[:, 2] - 0.5 * Lbox) ** 2
    > Rcloud ** 2
]
numPart_out = int(pos_out.shape[0])

# Aggregate all particles
pos = concatenate((pos_in, pos_out), axis=0)

numPart = int(pos.shape[0])

print(
    "The number of partciles in the cloud is %d, down from %d"
    % (numPart_in, numPart_in_side ** 3)
)
print(
    "The number of particles in the ambient medium is %d, down from %d"
    % (numPart_out, numPart_out_side ** 3)
)
print("The total number of particles is %d" % numPart)

h = ones(numPart) * 2.0 * Rcloud / numPart ** (1.0 / 3.0)

# Solid rotation for cloud particles
v = zeros((numPart, 3))
v[:numPart_in, 0] = -Omega * R * sin_phi
v[:numPart_in, 1] = Omega * R * cos_phi

# Other attributes
ids = linspace(1, numPart, numPart)
m = ones(numPart) * M / numPart_in

u = ones(numPart)
u[:numPart_in] *= P_in / ((gamma - 1) * rho_in)
u[numPart_in:] *= P_out / ((gamma - 1) * rho_out)

B = zeros((numPart, 3))
B[:, 2] = Bini

epsilon_lim = cbrt(M / (numPart_in * 1e-11)) / 3.086e18
print(epsilon_lim)
print(
    "The softening length you need to correctly resolve densities up to 1e-11 g cm^-3 is %f pc"
    % epsilon_lim
)

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
grp.create_dataset("Masses", data=m, dtype="f")
grp.create_dataset("SmoothingLength", data=h, dtype="f")
grp.create_dataset("InternalEnergy", data=u, dtype="f")
grp.create_dataset("ParticleIDs", data=ids, dtype="L")
grp.create_dataset("MagneticFluxDensities", data=B, dtype="f")

file.close()
