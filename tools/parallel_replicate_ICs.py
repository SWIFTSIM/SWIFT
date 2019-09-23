#!/usr/bin/env python
"""
Usage:
    python parallel_replicate_ICs.py IC_file.hdf5 rep_fac

where IC_file.hdf5 is the ICs file that you want to replicate and rep_fac is the
replication factor in each dimension

Description:
Reads in a ICs file and replicates the particles in each dimension by the
replication factor given and write a new IC called IC_file_xrep_fac.hdf5.

Example:
    python parallel_replicate_ICs.py EAGLE_ICs_50.hdf5 4

Running the above example will produce a tiled 50MPc box in each dimension to
give a 200MPc box.

Note: the script only replicates dark matter particles, but can be easily
extended to support other particle types.

This file is part of SWIFT.

Copyright (C) 2019 James Willis (james.s.willis@durham.ac.uk)
All Rights Reserved.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import h5py as h
import numpy as np
import matplotlib

matplotlib.use("Agg")
from pylab import *
import os.path
from tqdm import tqdm
from tqdm import trange
import time
from numba import jit, prange
from swiftsimio import Writer
from swiftsimio.units import cosmo_units

replicate = 1
box_size = 1
num_parts = 1

@jit(nopython=True, nogil=True)
def shift_pos(pos, pos_orig, i, j, k):
    
    offset = i * replicate * replicate + j * replicate + k

    # Copy original particle positions
    pos[offset * num_parts:(offset + 1) * num_parts] = pos_orig 

    # Shift positions
    shift = [i * box_size, j * box_size, k * box_size]

    for n in range(offset * num_parts, (offset + 1) * num_parts):
        pos[n][0] += shift[0]
        pos[n][1] += shift[1]
        pos[n][2] += shift[2]

@jit(nopython=True, parallel=True, nogil=True)
def parallel_replicate(pos, pos_orig):
    for i in prange(0,replicate):
        for j in prange(0,replicate):
            for k in prange(0,replicate):
                shift_pos(pos, pos_orig, i, j, k)
 
def main():

    # Parse command line arguments
    if len(sys.argv) < 3:
        print("Error: pass input file and replication factor (integer) as arguments.")
        print("python replicate_ICs.py EAGLE_ICs_50.hdf5 4")
        sys.exit()
    else:
        inputFile = sys.argv[1]
        global replicate
        replicate = int(sys.argv[2])
        if os.path.exists(inputFile) != 1:
            print("\n{} does not exist!\n".format(inputFile1))
            sys.exit()
  
    # Open ICs
    ics_file = h.File(inputFile, "r")

    replicate_factor = replicate * replicate * replicate

    global box_size, num_parts
    box_size = ics_file["/Header"].attrs["BoxSize"]
    num_parts = ics_file["/Header"].attrs["NumPart_Total"][1]

    print("Box size: {}".format(box_size))
    print("No. of original particles: {}".format(num_parts))
    print("New box size: {}".format(box_size * replicate))
    print("No. of replicated particles: {}".format(num_parts * replicate_factor))

    # Read input file fields
    pos_orig = ics_file["/PartType1/Coordinates"][:, :]
    mass_orig = ics_file["/PartType1/Masses"][:][0]
    
    # Create new arrays
    global pos, vel, mass, u, ids
    vel = pos = zeros((num_parts * replicate_factor, 3))
    ids = linspace(1, num_parts * replicate_factor, num_parts * replicate_factor)
    mass = zeros(num_parts * replicate_factor)
    mass[:] = mass_orig
    u = smoothing_length = zeros(num_parts * replicate_factor)

    start = time.time()
    # Replicate particles
    parallel_replicate(pos, pos_orig)

    print("Replicating particles took: %.3ss." % (time.time() - start))

    start = time.time()
    
    # Create output file
    base_filename = os.path.basename(inputFile)
    filename, file_extension = os.path.splitext(base_filename)
    outputFile = filename + "_x" + str(replicate) + ".hdf5"
    
    out_file = h.File(outputFile, 'w')

    # Copy Header and set new values
    ics_file.copy("/Header", out_file)
    grp = out_file["/Header"]
    grp.attrs["BoxSize"] = box_size * replicate
    grp.attrs["NumPart_Total"] =  [0, num_parts * replicate_factor, 0, 0, 0, 0]
    grp.attrs["NumPart_ThisFile"] = [0, num_parts * replicate_factor, 0, 0, 0, 0]
    
    # Copy Units
    ics_file.copy("/Units", out_file)

    # Particle group
    grp = out_file.create_group("/PartType1")
    grp.create_dataset('Coordinates', data=pos, dtype='d')
    grp.create_dataset('Velocities', data=vel, dtype='f')
    grp.create_dataset('Masses', data=mass, dtype='f')
    grp.create_dataset('SmoothingLength', data=smoothing_length, dtype='f')
    grp.create_dataset('InternalEnergy', data=u, dtype='f')
    grp.create_dataset('ParticleIDs', data=ids, dtype='L')
    
    print("Writing output file took: %.3ss." % (time.time() - start))

if __name__=="__main__":
    main()
