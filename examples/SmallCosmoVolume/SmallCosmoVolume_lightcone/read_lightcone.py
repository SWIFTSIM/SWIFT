#!/bin/env python

import h5py
import numpy as np
import healpy as hp


def read_map(basedir, basename, shell_nr, map_name, return_sum=False):
    """
    Read the specified healpix map for a lightcone shell
    """

    # Open the index file to determine number of files to read
    fname = "%s/%s_index.hdf5" % (basedir, basename)
    with h5py.File(fname, "r") as infile:
        nr_files_per_shell = infile["Lightcone"].attrs["nr_files_per_shell"][0]

    # Read the pixel data
    data = []
    for file_nr in range(nr_files_per_shell):
        fname = ("%s/%s_shells/shell_%d/%s.shell_%d.%d.hdf5" % 
                 (basedir, basename, shell_nr, basename, shell_nr, file_nr))
        with h5py.File(fname, "r") as infile:
            data.append(infile[map_name][...])
            if file_nr == 0 and return_sum:
                expected_sum = infile[map_name].attrs["expected_sum"]

    data = np.concatenate(data)

    if return_sum:
        return data, expected_sum
    else:
        return data


def read_particles(basedir, basename, part_type, properties):
    """
    Read particle data from a lightcone
    """

    # Open the index file to determine number of files to read
    fname = "%s/%s_index.hdf5" % (basedir, basename)
    with h5py.File(fname, "r") as infile:
        final_file_on_rank = infile["Lightcone"].attrs["final_particle_file_on_rank"]
        nr_mpi_ranks = infile["Lightcone"].attrs["nr_mpi_ranks"][0]

    # Make a dict to store the result
    data = {prop_name : [] for prop_name in properties}

    # Loop over MPI ranks
    for rank_nr in range(nr_mpi_ranks):
        # Loop over files written by this rank
        for file_nr in range(final_file_on_rank[rank_nr]+1):
            fname = "%s/%s_particles/%s_%04d.%d.hdf5" % (basedir, basename, basename, file_nr, rank_nr)
            with h5py.File(fname, "r") as infile:
                for prop_name in properties:
                    if part_type in infile:
                        data[prop_name].append(infile[part_type][prop_name][...])

    # Combine arrays from files
    for prop_name in properties:
        data[prop_name] = np.concatenate(data[prop_name])

    return data


