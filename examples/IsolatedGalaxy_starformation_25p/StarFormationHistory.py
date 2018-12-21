#!/usr/bin/env python3
import numpy as np
import h5py
import matplotlib.pyplot as plt

snaps = 114 

SFH = np.zeros(snaps)

for i in range(0,snaps):
    f = h5py.File('output_%04d.hdf5'%i, 'r')
    header = f["Header"]
    time = header.attrs['Time'][0]
    print(time)
    particles = f["PartType4"]
    newstar = particles["NewStarFlag"][:]
    SFH[i] = np.sum(newstar)

plt.plot(SFH)
plt.xlabel('Snapshot number')
plt.ylabel('Total number of stars formed')
plt.show()




f = h5py.File('output_%04d.hdf5'%i, 'r')
particles = f["PartType4"]
newstar = particles["NewStarFlag"][:]
Coordinates = particles["Coordinates"][:,:]
x = Coordinates[:,0][newstar==1]
y = Coordinates[:,1][newstar==1]

print(len(x), len(y))
plt.plot(x,y,'.')
plt.show()
