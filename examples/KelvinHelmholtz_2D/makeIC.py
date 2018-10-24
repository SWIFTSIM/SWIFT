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
import sys

# Generates a swift IC file for the Kelvin-Helmholtz vortex in a periodic box

# Parameters
L2    = 256       # Particles along one edge in the low-density region
gamma = 5./3.     # Gas adiabatic index
P1    = 2.5       # Central region pressure
P2    = 2.5       # Outskirts pressure
v1    = 0.5       # Central region velocity
v2    = -0.5      # Outskirts vlocity
rho1  = 2         # Central density
rho2  = 1         # Outskirts density
omega0 = 0.1
sigma = 0.05 / sqrt(2)
fileOutputName = "kelvinHelmholtz.hdf5" 
#---------------------------------------------------

# Start by generating grids of particles at the two densities
numPart2 = L2 * L2
L1 = int(sqrt(numPart2 / rho2 * rho1))
numPart1 = L1 * L1

#print "N2 =", numPart2, "N1 =", numPart1
#print "L2 =", L2, "L1 = ", L1
#print "rho2 =", rho2, "rho1 =", (float(L1*L1)) / (float(L2*L2))

coords1 = zeros((numPart1, 3))
coords2 = zeros((numPart2, 3))
h1 = ones(numPart1) * 1.2348 / L1
h2 = ones(numPart2) * 1.2348 / L2
m1 = zeros(numPart1)
m2 = zeros(numPart2)
u1 = zeros(numPart1)
u2 = zeros(numPart2)
vel1 = zeros((numPart1, 3))
vel2 = zeros((numPart2, 3))

# Particles in the central region
for i in range(L1):
    for j in range(L1):

        index = i * L1 + j
        
        x = i / float(L1) + 1. / (2. * L1)
        y = j / float(L1) + 1. / (2. * L1)

        coords1[index, 0] = x
        coords1[index, 1] = y
        u1[index] = P1 / (rho1 * (gamma-1.))
        vel1[index, 0] = v1
        
# Particles in the outskirts
for i in range(L2):
    for j in range(L2):

        index = i * L2 + j
        
        x = i / float(L2) + 1. / (2. * L2)
        y = j / float(L2) + 1. / (2. * L2)

        coords2[index, 0] = x
        coords2[index, 1] = y
        u2[index] = P2 / (rho2 * (gamma-1.))
        vel2[index, 0] = v2


# Now concatenate arrays
where1 = abs(coords1[:,1]-0.5) < 0.25
where2 = abs(coords2[:,1]-0.5) > 0.25

coords = append(coords1[where1, :], coords2[where2, :], axis=0)

#print L2*(L2/2), L1*(L1/2)
#print shape(coords), shape(coords1[where1,:]), shape(coords2[where2,:])
#print shape(coords), shape(logical_not(coords1[where1,:])), shape(logical_not(coords2[where2,:]))

vel = append(vel1[where1, :], vel2[where2, :], axis=0)
h = append(h1[where1], h2[where2], axis=0)
m = append(m1[where1], m2[where2], axis=0)
u = append(u1[where1], u2[where2], axis=0)
numPart = size(h)
ids = linspace(1, numPart, numPart)
m[:] = (0.5 * rho1 + 0.5 * rho2) / float(numPart)

# Velocity perturbation
vel[:,1] = omega0 * sin(4*pi*coords[:,0]) * (exp(-(coords[:,1]-0.25)**2 / (2 * sigma**2)) + exp(-(coords[:,1]-0.75)**2 / (2 * sigma**2)))
            
#File
fileOutput = h5py.File(fileOutputName, 'w')

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [1., 1., 0.1]
grp.attrs["NumPart_Total"] =  [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFileOutputsPerSnapshot"] = 1
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
ds[()] = coords
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


