###############################################################################
 # This file is part of SWIFT.
 # Copyright (c) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
import urllib
import os.path
import numpy as np
import sys

outputName = "bigCosmoVolume.hdf5" 


#--------------------------------------------------
if len(sys.argv) != 3:
    print "Invalid number of arguments. Need to provide down-sampling factor [0.,1.] and number of copies (integer)"
    print "Example: python makeIC.py 0.8 4"
    exit()

downsample = float(sys.argv[1])
n_copy = int(sys.argv[2])

if n_copy < 1:
    print "Number of copy must be >1"
    exit()

if downsample > 1. or downsample <= 0.:
    print "Down-sampling factor must be in [0,1]"
    exit()

#--------------------------------------------------

# Download the tile
if (not os.path.isfile("tile.hdf5")):
    print "Downloading initial tile..."
    urllib.urlretrieve ("http://virgodb.cosma.dur.ac.uk/swift-webstorage/ICs/tile.hdf5", "tile.hdf5")
    print "Done."
else:
    print "Tile already exists. No need to download..."

# Read in the tile
inputFile = h5py.File("tile.hdf5", 'r+')

grp = inputFile["/Header"]
boxSize = grp.attrs["BoxSize"]
numPart = grp.attrs["NumPart_Total"][0]

coords = inputFile["/PartType0/Coordinates"][:,:]
v = inputFile["/PartType0/Velocities"][:,:]
m = inputFile["/PartType0/Masses"][:]
h = inputFile["/PartType0/SmoothingLength"][:]
u = inputFile["/PartType0/InternalEnergy"][:]
ids = np.array(range(np.size(u)), dtype='L') + 1

# Downsample
print "Downsampling..."
indices = np.array(range(np.size(ids)))
np.random.shuffle(indices)

numPart *= downsample
indices = indices < numPart

coords = coords[indices,:]
v = v[indices,:]
m = m[indices]
h = h[indices] / 1.825742 # Correct from Gadget defintion of h to physical definition
u = u[indices]
ids = ids[indices]

numPart = np.size(ids)

# Now replicate the tile
if n_copy > 1:

    print "Tiling..."

    coords_tile = np.copy(coords)
    v_tile = np.copy(v)
    m_tile = np.copy(m)
    h_tile = np.copy(h)
    u_tile = np.copy(u)
    ids_tile = np.copy(ids)

    coords = np.zeros((0,3))
    v = np.zeros((0,3))
    m = np.zeros(0)
    h = np.zeros(0)
    u = np.zeros(0)
    ids = np.zeros(0, dtype='L')

    count = 0

    for i in range(n_copy):
        for j in range(n_copy):
            for k in range(n_copy):
                coords = np.append(coords, coords_tile + np.array([ i * boxSize[0], j * boxSize[1], k * boxSize[2] ]), axis=0)
                v = np.append(v, v_tile, axis=0)
                m = np.append(m, m_tile)
                h = np.append(h, h_tile)
                u = np.append(u, u_tile)
                
                ids = np.append(ids, ids_tile + count*numPart)
                count+=1


    numPart *= n_copy**3
    boxSize *= n_copy

# Copy the tile out
file = h5py.File(outputName, 'w')

# Header
grp = file.create_group("/Header")
grp.attrs["BoxSize"] = boxSize
grp.attrs["NumPart_Total"] =  [numPart, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [numPart, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = 0
grp.attrs["Dimension"] = 3

#Runtime parameters
grp = file.create_group("/RuntimePars")
grp.attrs["PeriodicBoundariesOn"] = 1

#Particle group
grp = file.create_group("/PartType0")
grp.create_dataset('Coordinates', data=coords, dtype='d')
grp.create_dataset('Velocities', data=v, dtype='f')
grp.create_dataset('Masses', data=m, dtype='f')
grp.create_dataset('SmoothingLength', data=h, dtype='f')
grp.create_dataset('InternalEnergy', data=u, dtype='f')
grp.create_dataset('ParticleIDs', data=ids, dtype='L')



