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

# Generates a swift IC file for the Square test in a periodic box

# Parameters
L = 64            # Number of particles on the side 
gamma = 5./3.     # Gas adiabatic index
rho0 = 4          # Gas central density
rho1 = 1          # Gas outskirt density
P0 = 2.5          # Gas central pressure
P1 = 2.5          # Gas central pressure
vx = 142.3        # Random velocity for all particles 
vy = -31.4
fileOutputName = "square.hdf5"
#---------------------------------------------------

vol = 1.

numPart_out = L * L
numPart_in = L * L * rho0 / rho1 / 4

L_out = int(sqrt(numPart_out))
L_in = int(sqrt(numPart_in))

pos_out = zeros((numPart_out, 3))
for i in range(L_out):
    for j in range(L_out):
        index = i * L_out + j
        pos_out[index, 0] =  i / (float(L_out)) + 1./(2. * L_out)
        pos_out[index, 1] =  j / (float(L_out)) + 1./(2. * L_out)
h_out = ones(numPart_out) * (1. / L_out) * 1.2348
m_out = ones(numPart_out) * vol  * rho1 / numPart_out
u_out = ones(numPart_out) * P1 / (rho1 * (gamma - 1.))

pos_in = zeros((numPart_in, 3))
for i in range(L_in):
    for j in range(L_in):
        index = i * L_in + j
        pos_in[index, 0] =  0.25 + i / float(2. * L_in) + 1./(2. * 2. * L_in)
        pos_in[index, 1] =  0.25 + j / float(2. * L_in) + 1./(2. * 2. * L_in)
h_in = ones(numPart_in) * (1. / L_in) * 1.2348
m_in = ones(numPart_in) * 0.25 * vol * rho0 / numPart_in
u_in = ones(numPart_in) * P0 / (rho0 * (gamma - 1.))

# Remove the central particles 
select_out = logical_or(logical_or(pos_out[:,0] < 0.25 , pos_out[:,0] > 0.75), logical_or(pos_out[:,1] < 0.25, pos_out[:,1] > 0.75))
pos_out = pos_out[select_out, :]
h_out = h_out[select_out]
u_out = u_out[select_out]
m_out = m_out[select_out]

# Add the central region
pos = append(pos_out, pos_in, axis=0)
h = append(h_out, h_in, axis=0)
u = append(u_out, u_in)
m = append(m_out, m_in)
numPart = size(h)
ids = linspace(1, numPart, numPart)
vel = zeros((numPart, 3))
vel[:,0] = vx
vel[:,1] = vy
        

#File
fileOutput = h5py.File(fileOutputName, 'w')

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [vol, vol, 0.2]
grp.attrs["NumPart_Total"] =  [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 2

#Units
grp = fileOutput.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.
grp.attrs["Unit mass in cgs (U_M)"] = 1.
grp.attrs["Unit time in cgs (U_t)"] = 1.
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.

#Particle group
grp = fileOutput.create_group("/PartType0")
ds = grp.create_dataset('Coordinates', (numPart, 3), 'd')
ds[()] = pos
ds = grp.create_dataset('Velocities', (numPart, 3), 'f')
ds[()] = vel
ds = grp.create_dataset('Masses', (numPart, 1), 'f')
ds[()] = m.reshape((numPart,1))
ds = grp.create_dataset('SmoothingLength', (numPart,1), 'f')
ds[()] = h.reshape((numPart,1))
ds = grp.create_dataset('InternalEnergy', (numPart,1), 'f')
ds[()] = u.reshape((numPart,1))
ds = grp.create_dataset('ParticleIDs', (numPart,1), 'L')
ds[()] = ids.reshape((numPart,1))

fileOutput.close()


