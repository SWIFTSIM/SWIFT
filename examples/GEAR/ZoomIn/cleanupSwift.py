#!/usr/bin/env python3

# ./translate_particles.py filename output_name
from h5py import File
import sys
from shutil import copyfile

NPartType = 6
filename = sys.argv[-2]
out = sys.argv[-1]

copyfile(filename, out)

f = File(out, "a")

# read cosmological parameters
a = f["Cosmology"].attrs["Scale-factor"][0]
h = f["Cosmology"].attrs["h"][0]

# add little h in the boxsize
f["Header"].attrs["BoxSize"] = f["Header"].attrs["BoxSize"][0] * h

# avoid the snapshot to be taken for SWIFT
del f["Header"].attrs["Code"]

# Update the fields (name + cosmo)
for i in range(NPartType):
    name = "PartType%i" % i
    if name not in f:
        continue

    if "SmoothingLengths" in f[name]:
        grp = f[name + "/SmoothingLengths"]
        grp[:] *= 1.823

    # fix issue due to the name of densities
    fields = [
        ("Coordinates", "Coordinates", h),
        ("Masses", "Masses", h),
        ("Velocities", "Velocities", 1. / a**0.5),
        ("Density", "Densities", 1. / h**2),
        ("Entropies", "Entropies", 1.),
        ("InternalEnergy", "InternalEnergies", 1. / a**2),
        ("SmoothingLength", "SmoothingLengths", h),
        ("Metals", "SmoothedElementAbundances", 1.)
    ]

    # create links
    for field in fields:
        if field[1] in f[name] and field[0] not in f[name]:
            f[name + "/" + field[0]] = f[name + "/" + field[1]]

    # apply unit transform
    for field in fields:
        field_name = name + "/" + field[0]
        if field[0] in f[name]:
            f[field_name][:] *= field[2]


# update cosmological parameters in order to be consistent with GEAR
cosmo = f["Cosmology"].attrs
head = f["Header"].attrs
head["Redshift"] = float(cosmo["Redshift"])
head["OmegaLambda"] = cosmo["Omega_lambda"]
head["Omega0"] = cosmo["Omega_m"]
head["HubbleParam"] = cosmo["h"][0]
head["Time"] = float(a)

f.close()
