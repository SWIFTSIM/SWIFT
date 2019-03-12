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
import sys
from numpy import *

# Generates a swift IC file containing a cartesian distribution of particles
# at a constant density and pressure in a cubic box

# Parameters
periodic= 1           # 1 For periodic box
boxSize = 1.
L = int(sys.argv[1])  # Number of particles along one axis
N = int(sys.argv[2])  # Write N particles at a time to avoid requiring a lot of RAM
rho = 2.              # Density
P = 1.                # Pressure
gamma = 5./3.         # Gas adiabatic index
eta = 1.2349      # 48 ngbs with cubic spline kernel
fileName = "uniformBox_%d.hdf5"%L

#---------------------------------------------------
numPart = L**3
mass = boxSize**3 * rho / numPart
internalEnergy = P / ((gamma - 1.)*rho)

#---------------------------------------------------

n_iterations = numPart / N
remainder = numPart % N

print "Writing", numPart, "in", n_iterations, "iterations of size", N, "and a remainder of", remainder

#---------------------------------------------------

#File
file = h5py.File(fileName, 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = boxSize
grp.attrs["NumPart_Total"] =  [numPart % (long(1)<<32), 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [numPart / (long(1)<<32), 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0

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


#Particle group
grp = file.create_group("/PartType0")

# First create the arrays in the file
ds_v = grp.create_dataset('Velocities', (numPart, 3), 'f', chunks=True, compression="gzip")
ds_m = grp.create_dataset('Masses', (numPart,1), 'f', chunks=True, compression="gzip")
ds_h = grp.create_dataset('SmoothingLength', (numPart,1), 'f', chunks=True, compression="gzip")
ds_u = grp.create_dataset('InternalEnergy', (numPart,1), 'f', chunks=True, compression="gzip")
ds_id = grp.create_dataset('ParticleIDs', (numPart, 1), 'L', chunks=True, compression="gzip")
ds_x = grp.create_dataset('Coordinates', (numPart, 3), 'd', chunks=True, compression="gzip")

# Now loop and create parts of the dataset
offset = 0
for n in range(n_iterations):
    v  = zeros((N, 3))
    ds_v[offset:offset+N,:] = v
    v = zeros(1)

    m = full((N, 1), mass)
    ds_m[offset:offset+N] = m
    m = zeros(1)

    h = full((N, 1), eta * boxSize / L)
    ds_h[offset:offset+N] = h
    h = zeros(1)

    u = full((N, 1), internalEnergy)
    ds_u[offset:offset+N] = u
    u = zeros(1)

    ids = linspace(offset, offset+N, N, endpoint=False).reshape((N,1))
    ds_id[offset:offset+N] = ids + 1
    x      = ids % L;
    y      = ((ids - x) / L) % L;
    z      = (ids - x - L * y) / L**2;
    ids    = zeros(1)
    coords = zeros((N, 3))
    coords[:,0] = z[:,0] * boxSize / L + boxSize / (2*L)
    coords[:,1] = y[:,0] * boxSize / L + boxSize / (2*L)
    coords[:,2] = x[:,0] * boxSize / L + boxSize / (2*L)
    ds_x[offset:offset+N,:] = coords
    coords  = zeros((1,3))

    offset += N
    print "Done", offset,"/", numPart, "(%.1f %%)"%(100*(float)(offset)/numPart)

# And now, the remainder
v  = zeros((remainder, 3))
ds_v[offset:offset+remainder,:] = v
v = zeros(1)

m = full((remainder, 1), mass)
ds_m[offset:offset+remainder] = m
m = zeros(1)

h = full((remainder, 1), eta * boxSize / L)
ds_h[offset:offset+remainder] = h
h = zeros(1)

u = full((remainder, 1), internalEnergy)
ds_u[offset:offset+remainder] = u
u = zeros(1)

ids = linspace(offset, offset+remainder, remainder, endpoint=False).reshape((remainder,1))
ds_id[offset:offset+remainder] = ids + 1
x      = ids % L;
y      = ((ids - x) / L) % L;
z      = (ids - x - L * y) / L**2;
coords = zeros((remainder, 3))
coords[:,0] = z[:,0] * boxSize / L + boxSize / (2*L)
coords[:,1] = y[:,0] * boxSize / L + boxSize / (2*L)
coords[:,2] = x[:,0] * boxSize / L + boxSize / (2*L)
ods_x[offset:offset+remainder,:] = coords

print "Done", offset+remainder,"/", numPart




file.close()
