#!/usr/bin/env python

"""
Usage:
  ./create_virtual_snapshot.py SNAPSHOT_NAME

where SNAPSHOT_NAME is one of multiple snapshots in a multi-file snapshot.
SNAPSHOT_NAME is assumed to be of the format PREFIX.COUNTER.hdf5, where COUNTER
runs from 0 to the number of files minus one.

This script will create a new file, PREFIX.hdf5, which contains the same metadata
as the original multi-file snapshot, but presents itself as a single snapshot file,
using HDF5's virtual dataset feature. For any other tool reading the file, it will
look as if the new virtual file is a single snapshot, except for the attribute
"Header:Virtual", which will be set to 1.

This file is part of SWIFT.

Copyright (C) Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

import numpy as np
import h5py
import argparse
import os
import re

# pre-compile some regular expressions
filename_re = re.compile("(\S+)\.([0-9]+)\.hdf5\Z")
parttype_re = re.compile("PartType([0-6]+)\Z")
partname_re = re.compile("(\S+)Particles\Z")

# match PartTypeX to a human friendly group name
partname_dict = {
    "PartType0": "GasParticles",
    "PartType1": "DMParticles",
    "PartType2": "DMBackgroundParticles",
    "PartType3": "SinkParticles",
    "PartType4": "StarsParticles",
    "PartType5": "BHParticles",
    "PartType6": "NeutrinoParticles",
}

# parse the single command-line argument
argparser = argparse.ArgumentParser(
    "Create a virtual snapshot for the given multi-file snapshot."
)
argparser.add_argument(
    "input",
    action="store",
    help="File name of one of the files in the multi-file snapshot.",
)
argparser.add_argument(
    "--force",
    action="store_true",
    help="Forcefully overwrite the virtual snapshot if it already exists.",
)
args = argparser.parse_args()

# parse the input file name
try:
    prefix, file_nr = filename_re.match(args.input).groups()
except:
    raise RuntimeError(
        f"Could not decompose input filename {args.input}."
        " Make sure it has the format PREFIX.FILE_NR.hdf5!"
    )
file_nr = int(file_nr)

# open the input file and get some useful information:
# - the number of files
# - the index of that particular file (sanity check)
# - the (total) number of particles
try:
    with h5py.File(args.input, "r") as handle:
        nfile = handle["/Header"].attrs["NumFilesPerSnapshot"][0]
        ifile = handle["/Header"].attrs["ThisFile"][0]
        particle_counts = handle["/Header"].attrs["NumPart_Total_HighWord"][:]
        particle_counts = particle_counts.astype(np.int64)
        particle_counts <<= 32
        particle_counts += handle["/Header"].attrs["NumPart_Total"][:]
except:
    raise RuntimeError(f"Cannot open file {args.input}!")

# perform a sanity check on the file index
if ifile != file_nr:
    raise RuntimeError(
        f"File index ({ifile}) does not match file name index ({file_nr})."
        " Something is wrong with these files!"
    )

# check that we are actually dealing with a multi-file snapshot
# (in principle, calling this script on a single file snapshot will fail when parsing the file name)
if nfile == 1:
    print(
        "You are running this script on a single file snapshot. No virtual file needs to be created!"
    )
    exit(0)

# compose the new virtual file name
output_filename = f"{prefix}.hdf5"

# check that the file does not exist yet
if os.path.exists(output_filename):
    print(f"Output file {output_filename} already exists!")
    if args.force:
        print("Forcefully overwriting the existing file, as requested.")
    else:
        print("Not overwriting it. Use --force to overwrite the existing file.")
        exit(0)

# copy all the groups (and datasets) that do not contain particles
# they do not require any changes, or only require minor changes
try:
    with h5py.File(args.input, "r") as ifile, h5py.File(output_filename, "w") as ofile:
        for group in ifile.keys():
            if parttype_re.match(group) is None and partname_re.match(group) is None:
                ifile.copy(group, ofile)
except:
    raise RuntimeError(
        "Something went wrong while trying to copy over non-particle groups from"
        f" {args.input} to {output_filename}!"
    )

# update the header of the virtual snapshot file
# the virtual file presents itself as an actual single file snapshot containing all the particles
# however, we set Header:Virtual to 1 to distinguish it from such a snapshot
try:
    with h5py.File(output_filename, "r+") as handle:
        handle["/Header"].attrs["NumFilesPerSnapshot"] = np.array([1], dtype=np.int32)
        handle["/Header"].attrs["ThisFile"] = np.array([0], dtype=np.int32)
        handle["/Header"].attrs["Virtual"] = np.array([1], dtype=np.int32)
        handle["/Header"].attrs["NumPart_ThisFile"] = particle_counts
except:
    raise RuntimeError(
        f"Could not update header properties of output file {output_filename}!"
    )

# now process the particle datasets
# first, we loop over all the files to obtain the data type and shape of each dataset
# we also count the number of elements in the dataset belonging to each sub-file
particle_groups = {}
for ifile in range(nfile):
    thisfile = f"{prefix}.{ifile}.hdf5"
    try:
        with h5py.File(thisfile, "r") as handle:
            for group in handle.keys():
                if parttype_re.match(group) is not None:
                    if not group in particle_groups:
                        particle_groups[group] = {}
                    for dset in handle[group].keys():
                        this_shape = handle[group][dset].shape
                        if not dset in particle_groups[group]:
                            if len(this_shape) == 2:
                                shape = this_shape[1]
                            else:
                                shape = 1
                            particle_groups[group][dset] = {
                                "dtype": handle[group][dset].dtype,
                                "shape": shape,
                                "sizes": np.zeros(nfile, dtype=np.int64),
                            }
                        particle_groups[group][dset]["sizes"][ifile] = this_shape[0]
    except:
        raise RuntimeError(
            f"Something went wrong while retrieving group information from {thisfile}!"
        )

# now we have all the information to create the new virtual datasets
# these present themselves as if they are a normal dataset containing values for all the particles
# however, the data is provided by a virtual link to the corresponding dataset in the original sub-file
try:
    with h5py.File(output_filename, "r+") as handle:
        for group in particle_groups:
            newgroup = handle.create_group(group)
            for dset in particle_groups[group]:
                path = f"{group}/{dset}"
                dtype = particle_groups[group][dset]["dtype"]
                shape = particle_groups[group][dset]["shape"]
                sizes = particle_groups[group][dset]["sizes"]
                offsets = np.cumsum(sizes)
                totsize = offsets[-1]
                offsets -= sizes
                if shape == 1:
                    real_shape = (totsize,)
                else:
                    real_shape = (totsize, shape)
                layout = h5py.VirtualLayout(shape=real_shape, dtype=dtype)
                for ifile in range(nfile):
                    if shape == 1:
                        this_shape = (sizes[ifile],)
                    else:
                        this_shape = (sizes[ifile], shape)
                    layout[
                        offsets[ifile] : offsets[ifile] + sizes[ifile]
                    ] = h5py.VirtualSource(
                        f"{prefix}.{ifile}.hdf5", path, shape=this_shape
                    )
                newgroup.create_virtual_dataset(dset, layout)
            # also create the soft link to the dataset, now that we are processing it anyway
            handle[partname_dict[group]] = h5py.SoftLink(group)
except:
    raise RuntimeError(
        f"Something went wrong while setting up virtual datasets in {output_filename}!"
    )

# copy over the attributes of the particle groups and datasets
try:
    with h5py.File(args.input, "r") as ifile, h5py.File(output_filename, "r+") as ofile:
        for group in ifile.keys():
            if parttype_re.match(group) is not None:
                for attr in ifile[group].attrs:
                    ofile[group].attrs[attr] = ifile[group].attrs[attr]
                for dset in ifile[group].keys():
                    path = f"{group}/{dset}"
                    for attr in ifile[path].attrs:
                        ofile[path].attrs[attr] = ifile[path].attrs[attr]
except:
    raise RuntimeError(
        "Something went wrong while copying dataset attributes from"
        f" {args.input} to {output_filename}!"
    )

# finally: update the cell metadata
# again, the virtual file presents itself as a single file snapshot
# this means that the Cells/Files elements are all 0
# it also means we have to update the offsets for all cells with the number of particles in cells
# belonging to files with a lower index
# we first gather these offsets (they are the same as the dataset offsets)
particle_offsets = {}
for group in particle_groups:
    for dset in particle_groups[group]:
        particle_offsets[group] = particle_groups[group][dset]["sizes"]
        break
    particle_offsets[group] = (
        np.cumsum(particle_offsets[group], dtype=np.int64) - particle_offsets[group]
    )

# now we actually update the cell metadata
try:
    with h5py.File(output_filename, "r+") as handle:
        for group in particle_offsets:
            files = handle[f"Cells/Files/{group}"][:]
            handle[f"Cells/OffsetsInFile/{group}"][:] += particle_offsets[group][files]
            handle[f"Cells/Files/{group}"][:] = 0
except:
    raise RuntimeError(
        f"Something went wrong while updating cell metadata for {output_filename}!"
    )

print(f"Finished writing {output_filename}.")
