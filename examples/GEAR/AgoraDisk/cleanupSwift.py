#!/usr/bin/env python3

# ./translate_particles.py filename output_name
from h5py import File
import sys
from shutil import copyfile

NPartType = 1
filename = sys.argv[-2]
out = sys.argv[-1]

copyfile(filename, out)

f = File(out)

name = "PartType0/ElementAbundance"
if name in f:
    del f[name]

for i in range(NPartType):
    name = "PartType%i" % i
    if name not in f:
        continue

    grp = f[name + "/SmoothingLength"]
    grp[:] *= 1.823

f.close()
