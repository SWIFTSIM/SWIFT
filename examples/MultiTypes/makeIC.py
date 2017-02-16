###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk),
 #                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
import sys
from numpy import *

# Generates a swift IC file containing a cartesian distribution of DM particles
# with a density of 1

# Parameters
periodic= 1           # 1 For periodic box
boxSize = 1.
Lgas = int(sys.argv[1])  # Number of particles along one axis
rhoGas = 2.              # Density
P = 1.                   # Pressure
gamma = 5./3.            # Gas adiabatic index
eta = 1.2349             # 48 ngbs with cubic spline kernel
rhoDM = 1.
Ldm = int(sys.argv[2])  # Number of particles along one axis

massStars = 0.1
Lstars = int(sys.argv[3])  # Number of particles along one axis

fileName = "multiTypes.hdf5"

#---------------------------------------------------
numGas = Lgas**3
massGas = boxSize**3 * rhoGas / numGas
internalEnergy = P / ((gamma - 1.)*rhoGas)

numDM = Ldm**3
massDM = boxSize**3 * rhoDM / numDM

numStars = Lstars**3
massStars = massDM * massStars


#--------------------------------------------------

#File
file = h5py.File(fileName, 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = boxSize
grp.attrs["NumPart_Total"] =  [numGas, numDM, 0, 0, numStars, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numGas, numDM, 0, 0, numStars, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, massDM, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 3

#Runtime parameters
grp = file.create_group("/RuntimePars")
grp.attrs["PeriodicBoundariesOn"] = periodic

#Units
grp = file.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = 1.
grp.attrs["Unit mass in cgs (U_M)"] = 1.
grp.attrs["Unit time in cgs (U_t)"] = 1.
grp.attrs["Unit current in cgs (U_I)"] = 1.
grp.attrs["Unit temperature in cgs (U_T)"] = 1.


# Gas Particle group
grp = file.create_group("/PartType0")

v  = zeros((numGas, 3))
ds = grp.create_dataset('Velocities', (numGas, 3), 'f')
ds[()] = v
v = zeros(1)

m = full((numGas, 1), massGas)
ds = grp.create_dataset('Masses', (numGas,1), 'f')
ds[()] = m
m = zeros(1)

h = full((numGas, 1), eta * boxSize / Lgas)
ds = grp.create_dataset('SmoothingLength', (numGas,1), 'f')
ds[()] = h
h = zeros(1)

u = full((numGas, 1), internalEnergy)
ds = grp.create_dataset('InternalEnergy', (numGas,1), 'f')
ds[()] = u
u = zeros(1)

ids = linspace(0, numGas, numGas, endpoint=False).reshape((numGas,1))
ds = grp.create_dataset('ParticleIDs', (numGas, 1), 'L')
ds[()] = ids + 1
x      = ids % Lgas;
y      = ((ids - x) / Lgas) % Lgas;
z      = (ids - x - Lgas * y) / Lgas**2;
coords = zeros((numGas, 3))
coords[:,0] = z[:,0] * boxSize / Lgas + boxSize / (2*Lgas)
coords[:,1] = y[:,0] * boxSize / Lgas + boxSize / (2*Lgas)
coords[:,2] = x[:,0] * boxSize / Lgas + boxSize / (2*Lgas)
ds = grp.create_dataset('Coordinates', (numGas, 3), 'd')
ds[()] = coords





# DM Particle group
grp = file.create_group("/PartType1")

v  = zeros((numDM, 3))
ds = grp.create_dataset('Velocities', (numDM, 3), 'f')
ds[()] = v
v = zeros(1)

m = full((numDM, 1), massDM)
ds = grp.create_dataset('Masses', (numDM,1), 'f')
ds[()] = m
m = zeros(1)

ids = linspace(0, numDM, numDM, endpoint=False).reshape((numDM,1))
ds = grp.create_dataset('ParticleIDs', (numDM, 1), 'L')
ds[()] = ids + Lgas**3 + 1
x      = ids % Ldm;
y      = ((ids - x) / Ldm) % Ldm;
z      = (ids - x - Ldm * y) / Ldm**2;
coords = zeros((numDM, 3))
coords[:,0] = z[:,0] * boxSize / Ldm + boxSize / (2*Ldm)
coords[:,1] = y[:,0] * boxSize / Ldm + boxSize / (2*Ldm)
coords[:,2] = x[:,0] * boxSize / Ldm + boxSize / (2*Ldm)
ds = grp.create_dataset('Coordinates', (numDM, 3), 'd')
ds[()] = coords



# Star Particle group
grp = file.create_group("/PartType4")

v  = zeros((numStars, 3))
ds = grp.create_dataset('Velocities', (numStars, 3), 'f')
ds[()] = v
v = zeros(1)

m = full((numStars, 1), massStars)
ds = grp.create_dataset('Masses', (numStars,1), 'f')
ds[()] = m
m = zeros(1)

ids = linspace(0, numStars, numStars, endpoint=False).reshape((numStars,1))
ds = grp.create_dataset('ParticleIDs', (numStars, 1), 'L')
ds[()] = ids + Lgas**3 + 1
x      = ids % Ldm;
y      = ((ids - x) / Ldm) % Ldm;
z      = (ids - x - Ldm * y) / Ldm**2;
coords = zeros((numStars, 3))
coords[:,0] = z[:,0] * boxSize / Ldm + boxSize / (2*Ldm)
coords[:,1] = y[:,0] * boxSize / Ldm + boxSize / (2*Ldm)
coords[:,2] = x[:,0] * boxSize / Ldm + boxSize / (2*Ldm)
ds = grp.create_dataset('Coordinates', (numStars, 3), 'd')
ds[()] = coords


file.close()
