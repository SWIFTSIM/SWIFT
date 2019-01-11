################################################################################
# This file is part of SWIFT.
# Copyright (c) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
from numpy import *

# Parameters
T_i = 1.e4               # Initial temperature of the gas (in K)
gamma = 5./3.            # Gas adiabatic index
numPart_gas_1D = 24      # Number of particles along each dimension
numPart_stars_1D = 1    # Number of particles along each dimension
fileName = "stellar_evolution.hdf5"


# Some units
Mpc_in_m = 3.08567758e22
Msol_in_kg = 1.98848e30
Gyr_in_s = 3.08567758e19
mH_in_kg = 1.6737236e-27

# Some constants
kB_in_SI = 1.38064852e-23
G_in_SI = 6.67408e-11

# Box extent and density
x_min = -0.5 * Mpc_in_m * 3.2e-3
x_max = 0.5 * Mpc_in_m * 3.2e-3
rho_0 = mH_in_kg * 1.e6

# SI system of units
unit_l_in_si = Mpc_in_m
unit_m_in_si = Msol_in_kg * 1.e10
unit_t_in_si = Gyr_in_s
unit_v_in_si = unit_l_in_si / unit_t_in_si
unit_u_in_si = unit_v_in_si**2

# Total number of particles
numPart_gas = numPart_gas_1D**3
numPart_stars = numPart_stars_1D**3

#---------------------------------------------------

# Set box size and interparticle distance
boxSize = x_max - x_min
delta_x = boxSize / numPart_gas_1D

# Set the particle mass
m_i = boxSize**3 * rho_0 / (numPart_gas + numPart_stars)
print(m_i/Msol_in_kg)

# Build the arrays
coords = zeros((numPart_gas, 3))
v = zeros((numPart_gas, 3))
ids = linspace(1, numPart_gas, numPart_gas)
m = zeros(numPart_gas)
h = zeros(numPart_gas)
u = zeros(numPart_gas)
star_coords = zeros((numPart_stars, 3))
star_v = zeros((numPart_stars, 3))
star_ids = linspace(1, numPart_stars, numPart_stars)
star_m = zeros(numPart_stars)
star_h = zeros(numPart_stars)
star_u = zeros(numPart_stars)

# Set the particles on the left
for i in range(numPart_gas_1D):
  for j in range(numPart_gas_1D):
    for k in range(numPart_gas_1D):
      index = i * numPart_gas_1D**2 + j * numPart_gas_1D + k
      #q = x_min + (i + 0.5) * delta_x
      coords[index,0] = random.sample()*boxSize + x_min
      coords[index,1] = random.sample()*boxSize + x_min
      coords[index,2] = random.sample()*boxSize + x_min
      u[index] = kB_in_SI * T_i / (gamma - 1.) / mH_in_kg
      h[index] = 1.2348 * delta_x
      m[index] = m_i
      v[index,0] = 0.
      v[index,1] = 0.
      v[index,2] = 0.

for i in range(numPart_stars_1D):
  for j in range(numPart_stars_1D):
    for k in range(numPart_stars_1D):
      index = i * numPart_stars_1D**2 + j * numPart_stars_1D + k
      #q = x_min + (i + 0.5) * delta_x
      if numPart_stars_1D > 1:
      	star_coords[index,0] = random.sample()*boxSize + x_min
      	star_coords[index,1] = random.sample()*boxSize + x_min
      	star_coords[index,2] = random.sample()*boxSize + x_min
      else:
      	star_coords[index,0] = 0.
      	star_coords[index,1] = 0.
      	star_coords[index,2] = 0.
      star_u[index] = kB_in_SI * T_i / (gamma - 1.) / mH_in_kg
      star_h[index] = 1.2348 * delta_x
      star_m[index] = m_i
      star_v[index,0] = 0.
      star_v[index,1] = 0.
      star_v[index,2] = 0.

# Unit conversion
coords /= unit_l_in_si
v /= unit_v_in_si
m /= unit_m_in_si
h /= unit_l_in_si
u /= unit_u_in_si
star_coords /= unit_l_in_si
star_v /= unit_v_in_si
star_m /= unit_m_in_si
star_h /= unit_l_in_si
star_u /= unit_u_in_si

boxSize /= unit_l_in_si

#File
file = h5py.File(fileName, 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = [boxSize, boxSize, boxSize]
grp.attrs["NumPart_Total"] =  [numPart_gas, 0, 0, 0, numPart_stars, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart_gas, 0, 0, 0, numPart_stars, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 3

#Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 100. * unit_l_in_si
grp.attrs["Unit mass in cgs (U_M)"] = 1000. * unit_m_in_si
grp.attrs["Unit time in cgs (U_t)"] = 1. * unit_t_in_si
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

#Particle group
grp = file.create_group("/PartType0")
grp.create_dataset('Coordinates', data=coords, dtype='d')
grp.create_dataset('Velocities', data=v, dtype='f')
grp.create_dataset('Masses', data=m, dtype='f')
grp.create_dataset('SmoothingLength', data=h, dtype='f')
grp.create_dataset('InternalEnergy', data=u, dtype='f')
grp.create_dataset('ParticleIDs', data=ids, dtype='L')

#Star group
grp = file.create_group("/PartType4")
grp.create_dataset('Coordinates', data=star_coords, dtype='d')
grp.create_dataset('Velocities', data=star_v, dtype='f')
grp.create_dataset('Masses', data=star_m, dtype='f')
grp.create_dataset('SmoothingLength', data=star_h, dtype='f')
grp.create_dataset('InternalEnergy', data=star_u, dtype='f')
grp.create_dataset('ParticleIDs', data=star_ids, dtype='L')

file.close()

#import pylab as pl
#from mpl_toolkits.mplot3d import Axes3D
#
#fig = pl.figure()
#ax = Axes3D(fig)
#ax.scatter(coords[:,0], coords[:,1], coords[:,2], c="r", alpha=0.01)
#ax.scatter(star_coords[0], star_coords[1], star_coords[2], c="k", alpha=1)
#pl.show()
