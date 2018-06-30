"""
###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2017
#
# Josh Borrow (joshua.borrow@durham.ac.uk)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# -----------------------------------------------------------------------------
#
# This program is a set of helper functions that write particle data to a
# GADGET .hdf5 file *handle*. It should also function as a piece of
# documentation on the required attributes for SWIFT to function when running
# in GADGET compatibility mode.
# 
# Example Usage:
#
# import write_gadget as wg
# import h5py as h5
$ import numpy as np
#
#
# N_PARTICLES = 1000
#
#
# with h5.File("test.hdf5", "w") as f:
#    wg.write_header(
#        f,
#        boxsize=100.,
#        flag_entropy=1,
#        np_total=[N_PARTICLES, 0, 0, 0, 0, 0],
#        np_total_hw=[0, 0, 0, 0, 0, 0]
#    )
#
#    wg.write_runtime_pars(
#        f,
#        periodic_boundary=1,
#    )
#
#    wg.write_units(
#        f,
#        current=1.,
#        length=1.,
#        mass=1.,
#        temperature=1.,
#        time=1.
#    )
#
#    wg.write_block(
#       f,
#       part_type=0,
#       pos=np.array([np.arange(N_PARTICLES)] * 3).T,
#       vel=np.array([np.zeros(N_PARTICLES)] * 3).T,
#       ids=np.arange(N_PARTICLES),
#       mass=np.ones(N_PARTICLES,
#       int_energy=np.zeros(N_PARTICLES),
#       smoothing=np.ones(NP_ARTICLES),
#
###############################################################################
"""


def write_header(f, boxsize, flag_entropy, np_total, np_total_hw, other=False):
    """ Writes the "Header" section of the hdf5 file. The parameters in this
        function that are required are the ones that are required for SWIFT
        to function; note that the MassTable is **ignored** and that all
        particle masses should be placed into the particle data arrays.

        @param: f | file handle
            - the file handle of the hdf5 object (use h5py.File(filename, "w")
              to open a file handle of the correct type).

        @param boxsize | float / list (2D / 3D)
            - the boxsize. If a float is given it is assumed that the box is
              the same size in all directions.

        @param flag_entropy | int (0/1)
            - sets Flag_Entropy_ICs. This is a historical variable for cross
              compatibility with Gadget-3

        @param np_total | list (6D)
            - the total number of particles required of each type.

                Type/Index | Symbolic Type Name
               ------------|--------------------
                    0      |       Gas
                    1      |       Halo
                    2      |       Disk
                    3      |       Bulge
                    4      |       Stars
                    5      |       Bndry

        @param np_total_hw | list (6D)
            - the number of high-word particles in the file.


        @param other | dictionary | optional
            - a dictionary with any other parameters that you wish to pass into
              the file header. They will be passed such that the key is the
              name of the attribute in the hdf5 file.

    """

    # We'll first build a dictionary to iterate through.

    default_attributes = {
        "BoxSize" : boxsize,
        "Flag_Entropy_ICs" : flag_entropy,
        "NumPart_Total" : np_total,
        "NumPart_Total_HighWord" : np_total_hw,
        "NumFilesPerSnapshot" : 1,  # required for backwards compatibility
        "NumPart_ThisFile" : np_total, # Also required for bw compatibility
    }

    if other:
        attributes = dict(default_attributes, **other)
    else:
        attributes = default_attributes

    header = f.create_group("Header")

    # For some reason there isn't a direct dictionary interface (at least one
    # that is documented, so we are stuck doing this loop...

    for name, value in attributes.items():
        header.attrs[name] = value

    return


def write_runtime_pars(f, periodic_boundary, other=False):
    """ Writes the "RuntimeParams" section in the hdf5 file. The parameters in
        this function that are required are also required for SWIFT to function.
        If you wish to pass extra arguments into the runtime parameters you
        may do that by providing a dictionary to other.

        @param: f | file handle
            - the file handle of the hdf5 object (use h5py.File(filename, "w")
              to open a file handle of the correct type).

        @param: periodic_boundary | int (0/1)
            - the condition for periodic boundary conditions -- they are 'on'
              if this variable is 1, and off if it is 0. Note that SWIFT
              currently requires periodic boundary conditions to run (as of
              August 2017).

        @param other | dictionary | optional
            - a dictionary with any other parameters that you wish to pass into
              the RuntimePars. They will be passed such that the key is the
              name of the attribute in the hdf5 file.
    """

    # First build the dictionary

    default_attributes = {
        "PeriodicBoundariesOn" : periodic_boundary,
    }

    if other:
        attributes = dict(default_attributes, **other)
    else:
        attributes = default_attributes

    runtime = f.create_group("RuntimePars")

    for name, value in attributes.items():
        runtime.attrs[name] = value

    return


def write_units(f, current, length, mass, temperature, time, other=False):
    """ Writes the "RuntimeParams" section in the hdf5 file. The parameters in
        this function that are required are also required for SWIFT to function.
        If you wish to pass extra arguments into the runtime parameters you
        may do that by providing a dictionary to other.

        @param: f | file handle
            - the file handle of the hdf5 object (use h5py.File(filename, "w")
              to open a file handle of the correct type).

        @param: current | float
            - the current conversion factor in cgs units.

        @param: length | float
            - the length conversion factor in cgs units.

        @param: mass | float
            - the mass conversion factor in cgs units.

        @param: temperature | float
            - the temperature conversion factor in cgs units.

        @param: time | float
            - the time conversion factor in cgs units.

        @param: other | dictionary | optional
            - a dictionary with any other parameters that you wish to pass into
              the units attributes. They will be passed such that the key is
              the name of the attribute in the hdf5 file.
    """

    # First build the dictionary

    default_attributes = {
        "Unit current in cgs (U_I)": current,
        "Unit length in cgs (U_L)": length,
        "Unit mass in cgs (U_M)": mass,
        "Unit temperature in cgs (U_T)": temperature,
        "Unit time in cgs (U_t)": time,
    }

    if other:
        attributes = dict(default_attributes, **other)
    else:
        attributes = default_attributes

    units = f.create_group("Units")

    for name, value in attributes.items():
        units.attrs[name] = value

    return


def write_block(f, part_type, pos, vel, ids, mass, int_energy, smoothing, other=False):
    """ Writes a given block of data to PartType{part_type}. The above
        required parameters are required for SWIFT to run.

        @param: f | file handle
            - the file handle of the hdf5 object.

        @param part_type | int (0-5):
            - the identifiying number of the particle type.

                Type/Index | Symbolic Type Name
               ------------|--------------------
                    0      |       Gas
                    1      |       Halo
                    2      |       Disk
                    3      |       Bulge
                    4      |       Stars
                    5      |       Bndry

        @param pos | numpy.array
            - the array of particle positions with shape (n_particles, 3).

        @param vel | numpy.array
            - the array of particle velocities with shape (n_particles, 3).

        @param ids | numpy.array
            - the ids of the particles. Please note that particle IDs in
              SWIFT must be strictly positive. Shape (n_particles, 1)

        @param mass | numpy.array
            - the masses of the particles. In SWIFT MassTable in the header
              is ignored and so particle masses must be provided here. Shape
              (n_particles, 1)

        @param int_energy | numpy.array
            - the internal energies of the particles. Shape (n_particles, 1)

        @param smoothing | numpy.array
            - the smoothing lenghts of the individual particles. Please note 
              that this cannot be supplied only in the parameterfile and must
              be provided on a particle-by-particle basis in SWIFT. Shape
              (n_particles, 1)

        @param: other | dictionary | optional
            - a dictionary with any other parameters that you wish to pass into
              the particle data. They will be passed such that the key is
              the name of the dataset in the hdf5 file.
    """
    
    # Build the dictionary

    default_data = {
        "Coordinates" : pos,
        "Velocities" : vel,
        "ParticleIDs" : ids,
        "Masses" : mass,
        "InternalEnergy" : int_energy,
        "SmoothingLength" : smoothing,
    }

    if other:
        data = dict(default_data, **other)
    else:
        data = default_data

    particles = f.create_group("PartType" + str(part_type))

    for name, value in data.items():
        particles.create_dataset(name, data=value)

    return

