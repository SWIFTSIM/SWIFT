################################################################################
# This file is part of SWIFT.
# Copyright (c) 2015 Bert Vandenbroucke (bert.vandenbroucke@ugent.be)
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
import random
import numpy as np
import sys

# Generates a swift IC file with some density gradients, to check the gradient
# reconstruction

# Parameters
gamma = 5./3.     # Gas adiabatic index
gridtype = "cartesian"
if len(sys.argv) > 1:
    gridtype = sys.argv[1]

# stretched cartesian box ######################################################
if gridtype == "stretched":
    fileName = "Gradients_stretched.hdf5"
    factor = 8
    boxSize = [ 1.0 , 1.0/factor , 1.0/factor ]
    L = 20
    nx1 = factor*L/2
    ny1 = L
    nz1 = L
    numfac = 2.
    nx2 = int(nx1/numfac)
    ny2 = int(ny1/numfac)
    nz2 = int(nz1/numfac)
    npart = nx1*ny1*nz1 + nx2*ny2*nz2
    vol = boxSize[0] * boxSize[1] * boxSize[2]
    partVol1 = 0.5*vol/(nx1*ny1*nz1)
    partVol2 = 0.5*vol/(nx2*ny2*nz2)

    coords = np.zeros((npart,3))
    h = np.zeros((npart,1))
    ids = np.zeros((npart,1), dtype='L')
    idx = 0
    dcell = 0.5/nx1
    for i in range(nx1):
        for j in range(ny1):
            for k in range(nz1):
                coords[idx,0] = (i+0.5)*dcell
                coords[idx,1] = (j+0.5)*dcell
                coords[idx,2] = (k+0.5)*dcell
                h[idx] = 0.56/nx1
                ids[idx] = idx
                idx += 1
    dcell = 0.5/nx2
    for i in range(nx2):
        for j in range(ny2):
            for k in range(nz2):
                coords[idx,0] = 0.5+(i+0.5)*dcell
                coords[idx,1] = (j+0.5)*dcell
                coords[idx,2] = (k+0.5)*dcell
                h[idx] = 0.56/nx2
                ids[idx] = idx
                idx += 1

# cartesian box ################################################################
if gridtype == "cartesian":
    fileName = "Gradients_cartesian.hdf5"
    boxSize = [ 1.0 , 1.0 , 1.0 ]
    nx = 20
    npart = nx**3
    partVol = 1./npart
    coords = np.zeros((npart,3))
    h = np.zeros((npart,1))
    ids = np.zeros((npart,1), dtype='L')
    idx = 0
    dcell = 1./nx
    for i in range(nx):
        for j in range(nx):
            for k in range(nx):
                coords[idx,0] = (i+0.5)*dcell
                coords[idx,1] = (j+0.5)*dcell
                coords[idx,2] = (k+0.5)*dcell
                h[idx] = 1.12/nx
                ids[idx] = idx
                idx += 1

# random box ###################################################################
if gridtype == "random":
    fileName = "Gradients_random.hdf5"
    boxSize = [ 1.0 , 1.0 , 1.0 ]
    glass = h5py.File("../Glass/glass_50000.hdf5", "r")
    coords = np.array(glass["/PartType0/Coordinates"])
    npart = len(coords)
    partVol = 1./npart
    h = np.zeros((npart,1))
    ids = np.zeros((npart,1), dtype='L')
    for i in range(npart):
        h[i] = 0.019
        ids[i] = i

v = np.zeros((npart,3))
m = np.zeros((npart,1))
rho = np.zeros((npart,1))
u = np.zeros((npart,1))

for i in range(npart):
    rhox = coords[i,0]
    if coords[i,0] < 0.75:
        rhox = 0.75
    if coords[i,0] < 0.25:
        rhox = 1.-coords[i,0]
    rhoy = 1.+boxSize[1]-coords[i,1]
    if coords[i,1] < 0.75*boxSize[1]:
        rhoy = 1. + 0.25*boxSize[1]
    if coords[i,1] < 0.25*boxSize[1]:
        rhoy = 1.+coords[i,1]
    rhoz = 1.
    rho[i] = rhox + rhoy + rhoz
    P = 1.
    u[i] = P / ((gamma-1.)*rho[i])
    if gridtype == "stretched":
        if coords[i,0] < 0.5:
            m[i] = rho[i] * partVol1
        else:
            m[i] = rho[i] * partVol2
    else:
        m[i] = rho[i] * partVol

#File
file = h5py.File(fileName, 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = boxSize
grp.attrs["NumPart_Total"] =  [npart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [npart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]

#Particle group
grp = file.create_group("/PartType0")
ds = grp.create_dataset('Coordinates', (npart, 3), 'd')
ds[()] = coords
ds = grp.create_dataset('Velocities', (npart, 3), 'f')
ds[()] = v
ds = grp.create_dataset('Masses', (npart,1), 'f')
ds[()] = m
ds = grp.create_dataset('Density', (npart,1), 'd')
ds[()] = rho
ds = grp.create_dataset('SmoothingLength', (npart,1), 'f')
ds[()] = h
ds = grp.create_dataset('InternalEnergy', (npart,1), 'd')
ds[()] = u
ds = grp.create_dataset('ParticleIDs', (npart,1), 'L')
ds[()] = ids[:]

file.close()
