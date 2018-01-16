#!/usr/bin/env python3
import h5py as h5
import sys

# First snapshot
snap_filename = "coolingBox_0000.hdf5"

# Read the initial state of the gas
f = h5.File(snap_filename,'r')
he = f["/PartType0/HeDensity"][:]

if (he != 5.).any():
    print("Test Failed")
    print("Received", he)

else:
    print("Test Succeed")
