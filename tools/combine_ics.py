#!/usr/bin/env python
"""
Usage:
    combine_ics.py input_file.0.hdf5 merged_file.hdf5 gzip_level

This file combines Gadget-2 type 2 (i.e. hdf5) initial condition files
into a single file that can be digested by SWIFT. 
This has mainly be tested for DM-only (parttype1) files but also works
smoothly for ICs including gas. The special case of a mass-table for
the DM particles is handled. No unit conversions are applied nor are
any scale-factors or h-factors changed.
The script applies some compression and checksum filters to the output
to save disk space. 
The last argument `gzip_level` is used to specify the level of compression
to apply to all the fields in the file. Use 0 to cancel all coompression.
The default value is `4`.

This file is part of SWIFT.
Copyright (C) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)

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

import sys
import h5py as h5
import numpy as np

# Store the compression level
gzip_level = 4
if len(sys.argv) > 3:
    gzip_level = sys.argv[3]

# First, we need to collect some information from the master file
main_file_name = str(sys.argv[1])[:-7]
print("Merging snapshots files with name", main_file_name)
master_file_name = main_file_name + ".0.hdf5"
print("Reading master information from", master_file_name)
master_file = h5.File(master_file_name, "r")
grp_header = master_file["/Header"]

num_files = grp_header.attrs["NumFilesPerSnapshot"]
tot_num_parts = grp_header.attrs["NumPart_Total"]
tot_num_parts_high_word = grp_header.attrs["NumPart_Total_HighWord"]
entropy_flag = grp_header.attrs["Flag_Entropy_ICs"]
box_size = grp_header.attrs["BoxSize"]
time = grp_header.attrs["Time"]

# Combine the low- and high-words
tot_num_parts = tot_num_parts.astype(np.int64)
for i in range(6):
    tot_num_parts[i] += (np.int64(tot_num_parts_high_word[i]) << 32)

# Some basic information
print("Reading", tot_num_parts, "particles from", num_files, "files.")

# Check whether there is a mass table
DM_mass = 0.0
mtable = grp_header.attrs.get("MassTable")
if mtable is not None:
    DM_mass = grp_header.attrs["MassTable"][1]
if DM_mass != 0.0:
    print("DM mass set to", DM_mass, "from the header mass table.")
else:
    print("Reading DM mass from the particles.")


# Create the empty file
output_file_name = sys.argv[2]
output_file = h5.File(output_file_name, "w-")


# Header
grp = output_file.create_group("/Header")
grp.attrs["NumFilesPerSnapshot"] = 1
grp.attrs["NumPart_Total"] = tot_num_parts
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = tot_num_parts
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["BoxSize"] = box_size
grp.attrs["Flag_Entropy_ICs"] = entropy_flag
grp.attrs["Time"] = time

# Create the particle groups
if tot_num_parts[0] > 0:
    grp0 = output_file.create_group("/PartType0")
if tot_num_parts[1] > 0:
    grp1 = output_file.create_group("/PartType1")
if tot_num_parts[4] > 0:
    grp4 = output_file.create_group("/PartType4")
if tot_num_parts[5] > 0:
    grp5 = output_file.create_group("/PartType5")


# Helper function to create the datasets we need
def create_set(grp, name, size, dim, dtype):
    if dim == 1:
        grp.create_dataset(
            name,
            (size,),
            dtype=dtype,
            chunks=True,
            compression="gzip",
            compression_opts=gzip_level,
            shuffle=True,
            fletcher32=True,
            maxshape=(size,),
        )
    else:
        grp.create_dataset(
            name,
            (size, dim),
            dtype=dtype,
            chunks=True,
            compression="gzip",
            compression_opts=gzip_level,
            shuffle=True,
            fletcher32=True,
            maxshape=(size, dim),
        )


# Create the required datasets
if tot_num_parts[0] > 0:
    create_set(grp0, "Coordinates", tot_num_parts[0], 3, "d")
    create_set(grp0, "Velocities", tot_num_parts[0], 3, "f")
    create_set(grp0, "Masses", tot_num_parts[0], 1, "f")
    create_set(grp0, "ParticleIDs", tot_num_parts[0], 1, "l")
    create_set(grp0, "InternalEnergy", tot_num_parts[0], 1, "f")
    create_set(grp0, "SmoothingLength", tot_num_parts[0], 1, "f")

if tot_num_parts[1] > 0:
    create_set(grp1, "Coordinates", tot_num_parts[1], 3, "d")
    create_set(grp1, "Velocities", tot_num_parts[1], 3, "f")
    create_set(grp1, "Masses", tot_num_parts[1], 1, "f")
    create_set(grp1, "ParticleIDs", tot_num_parts[1], 1, "l")

if tot_num_parts[4] > 0:
    create_set(grp4, "Coordinates", tot_num_parts[4], 3, "d")
    create_set(grp4, "Velocities", tot_num_parts[4], 3, "f")
    create_set(grp4, "Masses", tot_num_parts[4], 1, "f")
    create_set(grp4, "ParticleIDs", tot_num_parts[4], 1, "l")

if tot_num_parts[5] > 0:
    create_set(grp5, "Coordinates", tot_num_parts[5], 3, "d")
    create_set(grp5, "Velocities", tot_num_parts[5], 3, "f")
    create_set(grp5, "Masses", tot_num_parts[5], 1, "f")
    create_set(grp5, "ParticleIDs", tot_num_parts[5], 1, "l")

# Heavy-lifting ahead. Leave a last message.
print("Datasets created in output file")


# Special case of the non-zero mass table
if DM_mass != 0.0:
    masses = np.ones(tot_num_parts[1], dtype=np.float) * DM_mass
    grp1["Masses"][:] = masses


# Cumulative number of particles read/written
cumul_parts = [0, 0, 0, 0, 0, 0]

# Loop over all the files that are part of the snapshots
for f in range(num_files):

    file_name = main_file_name + "." + str(f) + ".hdf5"
    file = h5.File(file_name, "r")
    file_header = file["/Header"]
    num_parts = file_header.attrs["NumPart_ThisFile"]

    print(
        "Copying data from file",
        f,
        "/",
        num_files,
        ": num_parts = [",
        num_parts[0],
        num_parts[1],
        num_parts[4],
        num_parts[5],
        "]",
    )
    sys.stdout.flush()

    # Helper function to copy data
    def copy_grp(name_new, name_old, ptype):
        full_name_new = "/PartType" + str(ptype) + "/" + name_new
        full_name_old = "/PartType" + str(ptype) + "/" + name_old
        output_file[full_name_new][
            cumul_parts[ptype] : cumul_parts[ptype] + num_parts[ptype]
        ] = file[full_name_old]

    def copy_grp_same_name(name, ptype):
        copy_grp(name, name, ptype)

    if num_parts[0] > 0:
        copy_grp_same_name("Coordinates", 0)
        copy_grp_same_name("Velocities", 0)
        copy_grp_same_name("Masses", 0)
        copy_grp_same_name("ParticleIDs", 0)
        copy_grp_same_name("InternalEnergy", 0)
        copy_grp_same_name("SmoothingLength", 0)

    if num_parts[1] > 0:
        copy_grp_same_name("Coordinates", 1)
        copy_grp_same_name("Velocities", 1)
        copy_grp_same_name("ParticleIDs", 1)
        if DM_mass == 0.0:  # Do not overwrite values if there was a mass table
            copy_grp_same_name("Masses", 1)

    if num_parts[4] > 0:
        copy_grp_same_name("Coordinates", 4)
        copy_grp_same_name("Velocities", 4)
        copy_grp_same_name("Masses", 4)
        copy_grp_same_name("ParticleIDs", 4)

    if num_parts[5] > 0:
        copy_grp_same_name("Coordinates", 5)
        copy_grp_same_name("Velocities", 5)
        copy_grp_same_name("Masses", 5)
        copy_grp_same_name("ParticleIDs", 5)

    cumul_parts[0] += num_parts[0]
    cumul_parts[1] += num_parts[1]
    cumul_parts[4] += num_parts[4]
    cumul_parts[5] += num_parts[5]
    file.close()

print("All done! SWIFT is waiting.")
