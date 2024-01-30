#!/usr/bin/env python3

import h5py
import numpy as np
import unyt
from swiftsimio import Writer
from swiftsimio.units import cosmo_units

# Unit system we're working with
unitsystem = cosmo_units

# Box is 100 Mpc in each direction
boxsize = 260 * unyt.Mpc
boxsize = boxsize.to(unitsystem["length"])

reduced_speed_of_light_fraction = 1.

# Number of photon groups
nPhotonGroups = 1

# Number of particles in each dimension
n_p = 1000

# Filename of ICs to be generated
outputfilename = "uniform_3D.hdf5"

def initial_condition():
    """
    The initial conditions of the uniform box

    x: particle position. 3D unyt array
    unitsystem: Currently used unitsystem

    returns:
    E: Photon energy density
    F: Photon flux
    """

    uL = 3.0857E+24 # 1 Mpc in cm
    uT = 977.792221513146 * 365. * 24. * 3600. * 1e9 # *977.792221513146*Gyr in s
    uM = 1.98892e43 # 10000000000.0*Msun

    uE = uM * uL**2 / uT**2

    c_internal = 2.998e10 / (uL/uT) * reduced_speed_of_light_fraction

    # assume energies below are given in 1e10erg
    unit_conversion = 1e50 / uE
    
    # Uniform energy
    E = np.ones((n_p), dtype=np.float64) * unit_conversion * 1e5

    # Assuming all photons flow in only one direction
    # (optically thin regime, "free streaming limit"),
    # we have that |F| = c * E
    fluxes = np.zeros((3,n_p), dtype=np.float64)
    fluxes[0] *= E * c_internal / 1.73205  # sqrt(3)
    fluxes[1] *= E * c_internal / 1.73205  # sqrt(3)
    fluxes[2] *= E * c_internal / 1.73205  # sqrt(3)

    return E, fluxes.T

if __name__ in ("__main__"):
    coords = np.random.uniform(size = (n_p, 3)) * boxsize

    w = Writer(unitsystem, boxsize, dimension=3)

    w.gas.coordinates = coords
    w.gas.velocities = np.zeros((n_p, 3)) * (unyt.cm / unyt.s)

    mpart = 1e20 * unyt.M_sun
    mpart = mpart.to(unitsystem["mass"])
    w.gas.masses = np.ones(n_p, dtype=np.float64) * mpart
    w.gas.internal_energy = (np.ones(n_p, dtype=np.float64) * (300.0 * unyt.kb * unyt.K) / unyt.g)

    # Generate initial guess for smoothing lengths based on MIPS
    w.gas.generate_smoothing_lengths(boxsize=boxsize, dimension=3)

    # If IDs are not present, this automatically generates
    w.write(outputfilename)

    # Now open file back up again and add photon groups
    # you can make them unitless, the units have already been
    # written down in the writer, In this case, it's cosmo_units

    F = h5py.File(outputfilename, "r+")
    header = F["Header"]
    nparts = header.attrs["NumPart_ThisFile"][0]
    parts = F["/PartType0"]
    
    # Generate initial conditions
    E, fluxes = initial_condition()

    # Create photon energy data entry
    dsetname = "PhotonEnergiesGroup1"
    energydata = np.zeros((nparts), dtype=np.float32)
    parts.create_dataset(dsetname, data=E)
    
    # Create photon fluxes data entry
    dsetname = "PhotonFluxesGroup1"
    fluxdata = np.zeros((nparts, 3), dtype=np.float32)
    parts.create_dataset(dsetname, data=fluxes)


    
