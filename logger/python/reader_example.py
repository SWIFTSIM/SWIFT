#!/usr/bin/env python3
"""
Read a logger file by using an index.
Example: ./reader_example.py ../../examples/SedovBlast_3D/index.dump ../../examples/SedovBlast_3D/index_0005.hdf5
"""
import sys
from h5py import File
import numpy as np
import matplotlib.pyplot as plt
sys.path.append("../.libs/")

import liblogger as logger

# Get filenames
if len(sys.argv) != 3:
    print("WARNING missing arguments. Will use the default ones")
    index = "../../examples/SedovBlast_3D/index_0002.hdf5"
    dump = "../../examples/SedovBlast_3D/index.dump"
else:
    index = sys.argv[-1]
    dump = sys.argv[-2]

# constant
offset_name = "PartType0/Offset"
header = "Header"
time_name = "Time Offset"

# Read index file
with File(index, "r") as f:
    if offset_name not in f:
        raise Exception("Unable to find the offset dataset")
    offset = f[offset_name][:]

    if header not in f:
        raise Exception("Unable to find the header")
    if time_name not in f[header].attrs:
        raise Exception("Unable to find the time offset")
    time_offset = f[header].attrs[time_name]

# read dump
data = logger.loadFromIndex(offset, dump, time_offset)

# Compute distance from center
pos = data["positions"]
center = pos.mean()
r2 = np.sum((pos - center)**2, axis=1)

# plot entropy vs distance
plt.plot(np.sqrt(r2), data["entropy"], '.')

plt.xlim(0., 0.5)
plt.ylim(-5, 50)
plt.xlabel("Radius")
plt.ylabel("Entropy")
plt.show()
