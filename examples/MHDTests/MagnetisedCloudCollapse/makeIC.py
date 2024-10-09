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

# Constants
G   = 4.300918e-03
mu0 = 1.950089e-06
c1  = 0.53

# Parameters (taken from Hopkins 2016)
Rcloud = 4.628516371e16
Lbox = 10.0 * Rcloud

gamma = 4.0 / 3.0  # Gas adiabatic index

M = 1.99e33                          # total mass of the sphere
T = 4.7e5                            # initial orbital period in years
Omega = 2 * pi / (T * 3.1536e7)      # initial angular frequency of cloud 

mu   = 10                            # mass to magnetic field flux through sphere, normalised to a critical value for collapse. Refer to e.g. Henebelle & Fromang 2008 for details.
Bini = 3.0 / c1 * sqrt(mu0 * G / 5.0) * M / (Rcloud * Rcloud) * 1 / mu

# Barotropic EoS parameters
cs0 = 2e4
inv_rho_c = 1e14

# Attributes of cloud and ambient medium
volume_cloud = (4 / 3) * pi * Rcloud ** 3
volume_cloud_box = (2 * Rcloud) ** 3
volume_sim_box = Lbox ** 3

rho_in = M / volume_cloud
rho_out_to_rho_in = 1 / 360
rho_out = rho_out_to_rho_in * rho_in

P_in  = rho_in * cs0 * cs0 * sqrt(1.0 + (rho_in * inv_rho_c) ** gamma)
P_out = rho_out * cs0 * cs0 * sqrt(1.0 + (rho_out * inv_rho_c) ** gamma)

# Read glass files
fileName = "magnetised_cloud.hdf5"

glass = h5py.File("glassCube_128.hdf5", "r")
pos_gf = glass["/PartType0/Coordinates"][:, :]
h_gf   = glass["/PartType0/SmoothingLength"][:]

# Position cloud and ambient medium particles
cloud_box_side = 2.0 *	Rcloud
atmosphere_box_side = (1.0 / cbrt(rho_out_to_rho_in)) * cloud_box_side

pos_in = cloud_box_side * pos_gf 
h_in   = cloud_box_side * h_gf

pos_in -= 0.5 * cloud_box_side

mask_in = pos_in[:,0]**2 + pos_in[:,1]**2 + pos_in[:,2]**2 < Rcloud * Rcloud

pos_in = pos_in[mask_in]
h_in = h_in[mask_in]

numPart_in = int(len(h_in))

pos_out = atmosphere_box_side * pos_gf
h_out   = atmosphere_box_side * h_gf

pos_out -= 0.5 * atmosphere_box_side

mask_out = (pos_out[:,0]**2 + pos_out[:,1]**2 + pos_out[:,2]**2 > Rcloud * Rcloud) & (abs(pos_out[:,0]) < 0.5 * Lbox) & (abs(pos_out[:,1]) < 0.5 * Lbox) & (abs(pos_out[:,2]) < 0.5 * Lbox)

pos_out = pos_out[mask_out]
h_out = h_out[mask_out]

pos = concatenate((pos_in, pos_out), axis=0)
h = concatenate((h_in, h_out), axis=0)

numPart = int(len(h))

# Solid body rotation for cloud particles
mask = pos[:,0] ** 2 + pos[:,1] ** 2 + pos[:,2] ** 2 < Rcloud * Rcloud  

x = pos[:,0]
y = pos[:,1]

R = sqrt(x * x + y * y)

phi = arctan2(y, x)
cos_phi = cos(phi)
sin_phi = sin(phi)

v = zeros((numPart, 3))
v[mask][:, 0] = -Omega * R[mask] * sin_phi[mask]
v[mask][:, 1] = Omega * R[mask] * cos_phi[mask]

pos += 0.5 * Lbox

# Other attributes
ids = linspace(1, numPart, numPart)
m = ones(numPart) * M / numPart_in

u = ones(numPart)
u[mask] *= P_in / ((gamma - 1) * rho_in)
u[~mask] *= P_out / ((gamma - 1) * rho_out)

B = zeros((numPart, 3))
B[:, 2] = Bini

epsilon_lim = cbrt(M / (numPart_in * 1e-11)) / 3.086e18
print(epsilon_lim)
print(
    "The softening length you need to correctly resolve densities up to 1e-11 g cm^-3 is %f pc"
    % epsilon_lim
)

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
