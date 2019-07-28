"""
Script to convert the NIFTY ICs to those that are compatible with SWIFT.
Note that this leaves h-factors as-is to be fixed in-place by SWIFT.

You will need:

    + pygadgetreader from https://bitbucket.org/rthompson/pygadgetreader/overview
      (with 2to3 ran on all files to convert them to python3)
    + swiftsimio
"""

from pygadgetreader import *
from swiftsimio import Writer
import numpy as np
import unyt

filename = "/cosma7/data/dp004/jlvc76/nIFTy/IC_CLUSTER_00019"

length = unyt.kpc
mass = 1e10 * unyt.msun
time = unyt.s * unyt.kpc / unyt.km
velocity = length / time
energy_per_unit_mass = (length / time) ** 2


nifty_units = unyt.UnitSystem("nifty", length, mass, time)

writer = Writer(
    unit_system=nifty_units,
    box_size=readheader(filename, "boxsize") * length,
    dimension=3,
    compress=True,
    extra_header={
        "Redshift": readheader(filename, "redshift"),
        "Omega0": readheader(filename, "O0"),
        "OmegaLambda": readheader(filename, "Ol"),
        "HubbleParam": readheader(filename, "h"),
    },
)

writer.gas.coordinates = unyt.unyt_array(readsnap(filename, "pos", 0), length)

writer.gas.velocities = unyt.unyt_array(readsnap(filename, "vel", 0), velocity)

writer.gas.masses = unyt.unyt_array(readsnap(filename, "mass", 0), mass)

writer.gas.internal_energy = unyt.unyt_array(
    readsnap(filename, "u", 0), energy_per_unit_mass
)

# We must roll our own smoothing lengths.
n_part = len(writer.gas.masses)
x_range = writer.gas.coordinates.max() - writer.gas.coordinates.min()
mean_interparticle_sep = x_range / n_part ** (1 / 3)

writer.gas.smoothing_length = np.ones(n_part, dtype=float) * mean_interparticle_sep

writer.gas.particle_ids = unyt.unyt_array(readsnap(filename, "pid", 0), None)


def read_dm_quantity(name, unit, parttype):
    """
    The DM particles are in three sets because of their different masses.
    In SWIFT we have to combine these.
    """
    out = np.concatenate(
        [readsnap(filename, name, p) for p in parttype]
    ) * unit
    return out


for name, parttype in {"dark_matter": [1], "boundary": [2, 3, 5]}.items():
    writer_value = getattr(writer, name)

    writer_value.coordinates = read_dm_quantity("pos", length, parttype)

    writer_value.velocities = read_dm_quantity("vel", velocity, parttype)

    writer_value.masses = read_dm_quantity("mass", mass, parttype)

    writer_value.particle_ids = read_dm_quantity("pid", 1, parttype)

writer.write("nifty.hdf5")
