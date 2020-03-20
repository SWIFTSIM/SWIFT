#!/usr/bin/env

from h5py import File
from sys import argv

f = File(argv[-1], "r")

a = 0
# a = f["Cosmology"].attrs["Scale-factor"]
data = f["PartType0/Metals"][:]
print(a, data.min(), data.max())
