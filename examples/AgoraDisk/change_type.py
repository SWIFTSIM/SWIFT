#!/usr/bin/env python3

from h5py import File
from sys import argv
import numpy as np

"""
Change particle types in order to match the implemented types
"""

# number of particle type
N_type = 6


def getOption():
    if len(argv) != 2:
        raise IOError("You need to provide a filename")

    # get filename and read it
    filename = argv[-1]

    return filename


def groupName(part_type):
    return "PartType%i" % part_type


def changeType(f, old, new):
    # check if directory exists
    old_group = groupName(old)
    if old_group not in f:
        raise IOError("Cannot find group '%s'" % old)
    old = f[old_group]

    new = groupName(new)
    if new not in f:
        f.create_group(new)
    new = f[new]

    for name in old:
        tmp = old[name][:]
        del old[name]
        if name in new:
            tmp = np.append(new[name][:], tmp)
            del new[name]

        new.create_dataset(name, tmp.shape)
        new[name][:] = tmp

    del f[old_group]


def countPart(f):
    npart = []
    for i in range(N_type):
        name = groupName(i)
        if name in f:
            grp = f[groupName(i)]
            N = grp["Masses"].shape[0]
        else:
            N = 0
        npart.append(N)

    f["Header"].attrs["NumPart_ThisFile"] = npart
    f["Header"].attrs["NumPart_Total"] = npart
    f["Header"].attrs["NumPart_Total_HighWord"] = [0]*N_type


if __name__ == "__main__":
    filename = getOption()

    f = File(filename)

    changeType(f, 2, 1)
    changeType(f, 3, 1)
    changeType(f, 4, 1)

    countPart(f)

    f.close()
