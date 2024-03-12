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

reduced_speed_of_light_fraction = 1.0

# Number of photon groups
nPhotonGroups = 1

# Number of particles in each dimension
# Total number of particles is thus n_p^3
n_p = 10

# Filename of ICs to be generated
outputfilename = "uniform_3D.hdf5"


def initial_condition(unitsystem):
    """
    The initial conditions of the uniform box
    
    unitsystem: The unit system to use for IC

    returns:
    E: Photon energy density
    F: Photon flux
    """
    # you can make the photon quantities unitless, the units will
    # already have been written down in the writer.
    # However, that means that you need to convert them manually.

    unit_energy = (
        unitsystem["mass"] * unitsystem["length"] ** 2 / unitsystem["time"] ** 2
    )
    unit_velocity = unitsystem["length"] / unitsystem["time"]
    unit_flux = unit_energy * unit_velocity

    c_internal = (unyt.c * reduced_speed_of_light_fraction).to(unit_velocity)

    # Uniform energy
    E = np.ones((n_p ** 3), dtype=np.float64) * unit_energy

    # Assuming all photons flow in only one direction
    # (optically thin regime, "free streaming limit"),
    # we have that |F| = c * E
    fluxes = np.zeros((3, n_p ** 3), dtype=np.float64)
    fluxes[0] *= (E * c_internal / 1.73205).to(unit_flux)  # sqrt(3)
    fluxes[1] *= (E * c_internal / 1.73205).to(unit_flux)  # sqrt(3)
    fluxes[2] *= (E * c_internal / 1.73205).to(unit_flux)  # sqrt(3)

    return E, fluxes.T


if __name__ in ("__main__"):
    # Coordinate array
    coords = np.zeros((n_p ** 3, 3), dtype=np.float64)

    # Calculate grid of evenly spaced coordinates
    coords_per_dim = np.linspace(0.5, n_p - 0.5, n_p)
    grid = np.meshgrid(coords_per_dim, coords_per_dim, coords_per_dim)

    for i in range(3):
        coords[:, i] = grid[i].flatten()

    # Calculate and apply grid spacing
    dx = boxsize / n_p
    coords *= dx

    w = Writer(unitsystem, boxsize, dimension=3)

    w.gas.coordinates = coords
    w.gas.velocities = np.zeros((n_p ** 3, 3)) * (unyt.cm / unyt.s)

    mpart = 1e20 * unyt.M_sun
    mpart = mpart.to(unitsystem["mass"])
    w.gas.masses = np.ones(n_p ** 3, dtype=np.float64) * mpart
    w.gas.internal_energy = (
        np.ones(n_p ** 3, dtype=np.float64) * (300.0 * unyt.kb * unyt.K) / unyt.g
    )

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
    E, fluxes = initial_condition(unitsystem)

    # Create photon energy data entry
    dsetname = "PhotonEnergiesGroup1"
    energydata = np.zeros((nparts), dtype=np.float32)
    parts.create_dataset(dsetname, data=E)

    # Create photon fluxes data entry
    dsetname = "PhotonFluxesGroup1"
    fluxdata = np.zeros((nparts, 3), dtype=np.float32)
    parts.create_dataset(dsetname, data=fluxes)
