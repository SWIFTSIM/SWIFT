#!/usr/bin/env python

# ./translate_particles.py filename
from h5py import File
import sys

NPartType = 1
f = File(sys.argv[-1])

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
