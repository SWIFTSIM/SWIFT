###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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

import h5py
from numpy import *

# Some constants
solar_mass_cgs = 1.988480e33 
kpc_in_cm = 3.085678e21
mp_cgs = 1.67e-24
boltzmann_k_cgs = 1.38e-16

# Parameters
gamma = 5./3.      				# Gas adiabatic index
rho_cgs = mp_cgs        			# Background density
u0_cgs = 1.2e12					# Desired initial internal energy (1.2e12 ~ 10^4K)
P_cgs = rho_cgs*u0_cgs*(gamma - 1.)          	# Background pressure
fileName = "bh_swallowing.hdf5" 

# Units
unit_l_cgs = 3.085678e24  # kpc
unit_m_cgs = 1.988480e43  # 10^10 Msun
unit_v_cgs = 1e5          # km / s
unit_A_cgs = 1.
unit_T_cgs = 1.
unit_t_cgs = unit_l_cgs / unit_v_cgs

boxsize_cgs = 10. * kpc_in_cm
vol_cgs = boxsize_cgs**3

#---------------------------------------------------
glass = h5py.File("glassCube_32.hdf5", "r")

# Read particle positions and h from the glass
pos = glass["/PartType0/Coordinates"][:,:]
eps = 1e-6
pos = (pos - pos.min()) / (pos.max() - pos.min() + eps) * boxsize_cgs / unit_l_cgs
h = glass["/PartType0/SmoothingLength"][:] * boxsize_cgs / unit_l_cgs

numPart = size(h)

# Generate extra arrays
v = zeros((numPart, 3))
ids = linspace(1, numPart, numPart)
m_cgs = zeros(numPart)
u_cgs = zeros(numPart)
m = zeros(numPart)
u = zeros(numPart)

m_cgs[:] = rho_cgs * vol_cgs / numPart    
u_cgs[:] = P_cgs / (rho_cgs * (gamma - 1))

# BHs
bh_pos = zeros((2, 3))
bh_pos[0,:] = 0.5 * boxsize_cgs / unit_l_cgs

#diff = [4.812127e-03, 4.908179e-03, -4.878537e-03] - bh_pos[0,:]

#print diff

bh_pos[1,:] = 0.5 * boxsize_cgs / unit_l_cgs + diff

print bh_pos

# BHs don't move
bh_v = zeros((2, 3))
bh_v[:,:] = 0.

# increase mass to keep it at center
bh_m_cgs = ones(2) * m_cgs[0]
bh_ids = linspace(numPart + 1, numPart + 2, 2, dtype='L')
bh_h = ones(2) * [h.max()]

print bh_ids

#--------------------------------------------------

# Check quantities are correct for debugging
print("boxsize kpc " + str(boxsize_cgs/kpc_in_cm))
print("density cm^-3 " + str(rho_cgs/mp_cgs))
print("initial temperature K " + str(u_cgs[0] / boltzmann_k_cgs*((gamma - 1)*rho_cgs)))

# Convert to internal units
bh_m = bh_m_cgs/unit_m_cgs
m[:] = m_cgs/unit_m_cgs
u[:] = u_cgs*unit_v_cgs**-2
boxsize = boxsize_cgs/unit_l_cgs


#--------------------------------------------------

#File
file = h5py.File(fileName, 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = [boxsize]*3
grp.attrs["NumPart_Total"] =  [numPart, 0, 0, 0, 0, 2]
#grp.attrs["NumPart_Total"] =  [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 2]
#grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 3

#Runtime parameters
grp = file.create_group("/RuntimePars")
grp.attrs["PeriodicBoundariesOn"] = 0

#Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = unit_l_cgs
grp.attrs["Unit mass in cgs (U_M)"] = unit_m_cgs
grp.attrs["Unit time in cgs (U_t)"] = unit_t_cgs
grp.attrs["Unit current in cgs (U_I)"] = unit_A_cgs
grp.attrs["Unit temperature in cgs (U_T)"] = unit_T_cgs

#Particle group
grp = file.create_group("/PartType0")
grp.create_dataset('Coordinates', data=pos, dtype='d')
grp.create_dataset('Velocities', data=v, dtype='f')
grp.create_dataset('Masses', data=m, dtype='f')
grp.create_dataset('SmoothingLength', data=h, dtype='f')
grp.create_dataset('InternalEnergy', data=u, dtype='f')
grp.create_dataset('ParticleIDs', data=ids, dtype='L')

# stellar group
grp = file.create_group("/PartType5")
grp.create_dataset("Coordinates", (2,3),data=bh_pos, dtype="d")
grp.create_dataset('Velocities', (2,3),data=bh_v, dtype='f')
grp.create_dataset('Masses', (2,1),data=bh_m, dtype='f')
grp.create_dataset('SmoothingLength', (2,1),data=bh_h, dtype='f')
grp.create_dataset('ParticleIDs', (2,1),data=bh_ids[::-1], dtype='L')

file.close()
