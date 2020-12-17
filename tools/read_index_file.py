#!/usr/bin/env python3

import sys
import numpy as np

filename = sys.argv[-1]
n_type = 6


# dtype for the particle's data
dt = np.dtype([("ids", np.ulonglong),
               ("offset", np.uint64)])

# Read the file
with open(filename, "rb") as f:
    # read the time
    time = np.fromfile(f, dtype=float, count=1)
    time_int = np.fromfile(f, dtype=np.longlong, count=1)
    print("Time: {}, integer time: {}".format(
        time[0], time_int[0]))

    # read the number of particles
    nparts = np.fromfile(f, dtype=np.uint64, count=n_type)

    print("Number of particles:", nparts)

    # read if the file is sorted
    sort = np.fromfile(f, dtype=np.bool, count=1)
    print("File is sorted?", sort[0])

    # read the memory alignment garbage
    n = ((f.tell() + 7) & ~7) - f.tell()
    f.read(n)

    # read the particles
    print("Particles data (ids / offset):")
    for n in nparts:
        if n == 0:
            continue

        data = np.fromfile(f, dtype=dt, count=n)

        print("\t", data)

    # print the history of new particles
    n_new = np.fromfile(f, dtype=np.uint64, count=n_type)
    print("New particles: ", n_new)

    for n in n_new:
        if n == 0:
            continue

        data = np.fromfile(f, dtype=dt, count=n)
        print("\t", data)

    # print the history of particles removed
    n_rem = np.fromfile(f, dtype=np.uint64, count=n_type)
    print("Particles removed: ", n_rem)

    for n in n_rem:
        if n == 0:
            continue

        data = np.fromfile(f, dtype=dt, count=n)
        print("\t", data)
