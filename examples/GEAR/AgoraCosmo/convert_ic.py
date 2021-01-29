#!/usr/env/python3

from h5py import File
import shutil
import numpy as np
import sys
from swiftsimio.visualisation.smoothing_length_generation import generate_smoothing_lengths
import time
from unyt import kpc

filename = "agora_arepo.hdf5"
output = "agora_swift.hdf5"

shutil.copy(filename, output)

with File(output, "a") as f:

    f["Header"].attrs["Flag_Entropy_ICs"] = 0

    npart = f["Header"].attrs["NumPart_ThisFile"]
    mass = f["Header"].attrs["MassTable"]
    box = f["Header"].attrs["BoxSize"]
    T = f["Header"].attrs["suggested_gas_Tinit"]

    # Create the units
    f.create_group("Units")
    print("Assuming defaults AREPO units")
    u_l = 3.085678e21
    u_m = 1.989e43
    u_v = 1e5
    u_t = u_l / u_v
    f["Units"].attrs["Unit length in cgs (U_L)"] = u_l
    f["Units"].attrs["Unit mass in cgs (U_M)"] = u_m
    f["Units"].attrs["Unit time in cgs (U_t)"] = u_t
    f["Units"].attrs["Unit current in cgs (U_I)"] = 1.
    f["Units"].attrs["Unit temperature in cgs (U_T)"] = 1.

    # Create the mass arrays
    for i in range(len(npart)):
        if npart[i] == 0:
            continue

        grp = f["PartType%i" % i]
        if "Masses" in grp:
            continue
        masses = np.ones(npart[i]) * mass[i]
        grp.create_dataset('Masses', data=masses, dtype='f')

    # Create the smoothing lengths
    pos = f["PartType0/Coordinates"][:] * kpc
    h = generate_smoothing_lengths(pos, box * kpc, kernel_gamma=1.825742)
    f["PartType0"].create_dataset("SmoothingLength", data=h.value, dtype="f")

    # Deal with the internal energy
    u = np.ones(h.shape) * -1  # We set it through the parameters => fill it with garbage
    f["PartType0"].create_dataset("InternalEnergy", data=u, dtype="f")
