#!/usr/bin/env python3
"""
Read a logger file by using an index file and then write a snapshot.
Example: python3 create_snapshot -t 0.1 -o out.hdf5 ../../examples/SedovBlast_3D/index_*dump
"""
import sys
import h5py
import numpy as np
from glob import glob
import argparse
sys.path.append("../.libs/")

import liblogger as logger


def parse_arguments():
    parser = argparse.ArgumentParser(
        description='Read a logfile and plots some basic properties')

    default_files = "../../examples/HydroTests/SedovBlast_3D/index_*dump"
    default_files = glob(default_files)

    parser.add_argument("-o", '--output', dest='out',
                        type=str, default="out.hdf5",
                        help='Output filename')
    parser.add_argument("-t", '--time', dest='time',
                        type=float, default=0.01,
                        help='Simulation time')
    parser.add_argument('files', metavar='filenames', type=str, nargs="*",
                        help='The filenames of the logfiles')
    args = parser.parse_args()
    if len(args.files) == 0:
        args.files = default_files
    return args


def write_particle_type(snap, part_type, args):
    """
    Write a particle type inside an HDF5 file.

    Parameters
    ----------

    snap: h5py.File
        The output file.

    part_type: int
        The particle type to write

    args: ArgumentParser
        The argument parser

    Returns
    -------

    npart: int
        The number of particles of the given type.
    """
    fields = None
    fields_name = None

    # Read all the fields
    for f in args.files:
        if f.endswith(".dump"):
            filename = f[:-5]
        else:
            raise Exception("It seems that you are not providing a logfile (.dump)")

        # Open the logger
        with logger.Reader(filename, verbose=0) as reader:

            # Get the list of fields
            if fields_name is None:
                fields_name = reader.get_list_fields(part_type)

            # Abort if only requesting a type not implemented
            if "SpecialFlags" in fields_name:
                print("Part type %i not implemented, skipping it" % part_type)
                continue

            # Read the fields
            fields_tmp = reader.get_particle_data(fields_name, args.time, part_type)
            if fields is None:
                fields = fields_tmp
            else:
                for i, field in enumerate(fields):
                    fields[i] = np.append(fields[i], fields_tmp[i], axis=0)

    # Do we have this particle type?
    if fields is None or fields[0].shape[0] == 0:
        return 0

    # Get the number of particles
    npart = fields[0].shape[0]

    # Create the group
    name = "PartType%i" % part_type
    grp = snap.create_group(name)

    # Save the data
    for i, field in enumerate(fields_name):
        grp.create_dataset(field, data=fields[i])

    return npart


if __name__ == "__main__":
    args = parse_arguments()
    print("Output: %s" % args.out)
    print("basename: %s" % args.files)
    print("time: %g" % args.time)

    # read the logger
    n_types = 6

    # Create a snapshot
    with h5py.File(args.out, "w") as snap:
        npart = np.zeros(n_types)
        for part_type in range(n_types):
            npart[part_type] = write_particle_type(snap, part_type, args)


        # Write the header
        grp = snap.create_group("Header")
        grp.attrs["NumPart_Total"] = npart
        grp.attrs["NumPart_Total_HighWord"] = [0] * n_types
        grp.attrs["NumPart_ThisFile"] = npart
        grp.attrs["Time"] = args.time
        grp.attrs["NumFilesPerSnapshot"] = 1
        grp.attrs["MassTable"] = [0.0] * n_types

