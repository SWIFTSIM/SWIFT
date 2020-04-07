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
    name = "PartType{}/ElementAbundances".format(i)
    if name in f:
        del f[name]

for i in range(NPartType):
    name = "PartType%i" % i
    if name not in f:
        continue

    grp = f[name + "/SmoothingLengths"]
    grp[:] *= 1.823

    # fix issue due to the name of densities
    fields = [("Density", "Densities"),
              ("Entropies", "Entropies"),
              ("InternalEnergy", "InternalEnergies"),
              ("SmoothingLength", "SmoothingLengths")]

    for field in fields:
        if field[1] in f[name] and field[0] not in f[name]:
            f[name + "/" + field[0]] = f[name + "/" + field[1]]




cosmo = f["Cosmology"].attrs
head = f["Header"].attrs
head["OmegaLambda"] = cosmo["Omega_lambda"]
head["Omega0"] = cosmo["Omega_b"]
head["HubbleParam"] = cosmo["H0 [internal units]"]
head["Time"] = head["Time"][0]

f.close()
