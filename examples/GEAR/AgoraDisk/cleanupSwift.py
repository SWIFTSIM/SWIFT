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

for i in range(6):
    name = "PartType{}/ElementAbundance".format(i)
    if name in f:
        del f[name]

for i in range(NPartType):
    name = "PartType%i" % i
    if name not in f:
        continue

    grp = f[name + "/SmoothingLength"]
    grp[:] *= 1.823

cosmo = f["Cosmology"].attrs
head = f["Header"].attrs
head["OmegaLambda"] = cosmo["Omega_lambda"]
head["Omega0"] = cosmo["Omega_b"]
head["HubbleParam"] = cosmo["H0 [internal units]"]

f.close()
