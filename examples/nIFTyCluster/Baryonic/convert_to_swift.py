"""
Script to convert the NIFTY ICs to those that are compatible with SWIFT.
Note that this leaves h-factors as-is to be fixed in-place by SWIFT.

You will need:

    + pygadgetreader from https://bitbucket.org/rthompson/pygadgetreader/overview
      (with 2to3 ran on all files to convert them to python3)
    + swiftsimio
"""

from pygadgetreader import *
import swiftsimio as sw
import numpy as np
import unyt

filename = "/cosma7/data/dp004/jlvc76/nIFTy/IC_CLUSTER_00019"

length = unyt.kpc
mass = 1e10 * unyt.msun
time = (1.0 * unyt.s * unyt.kpc / unyt.km).to("s")
velocity = length / time
energy_per_unit_mass = (length / time) ** 2


nifty_units = unyt.UnitSystem("nifty", 1e3 * length, mass, time)

bs = readheader(filename, "boxsize") * length
boxsize_cosmo = sw.cosmo_array(
    [bs.value, bs.value, bs.value],
    bs.units,
    comoving=True,
    scale_factor=1.0,
    scale_exponent=1,
)
writer = sw.Writer(
    unit_system=nifty_units,
    boxsize=boxsize_cosmo,
    dimension=3,
    compress=True,
    extra_header={
        "Redshift": readheader(filename, "redshift"),
        "Omega0": readheader(filename, "O0"),
        "OmegaLambda": readheader(filename, "Ol"),
        "HubbleParam": readheader(filename, "h"),
    },
)

gas_coords = unyt.unyt_array(readsnap(filename, "pos", 0), length)
writer.gas.coordinates = sw.cosmo_array(
    gas_coords.value,
    gas_coords.units,
    comoving=True,
    scale_factor=1.0,
    scale_exponent=1,
)

gas_vel = unyt.unyt_array(readsnap(filename, "vel", 0), velocity)
writer.gas.velocities = sw.cosmo_array(
    gas_vel.value, gas_vel.units, comoving=True, scale_factor=1.0, scale_exponent=0
)

gas_mass = unyt.unyt_array(readsnap(filename, "mass", 0), mass)
writer.gas.masses = sw.cosmo_array(
    gas_mass.value, gas_mass.units, comoving=True, scale_factor=1.0, scale_exponent=0
)

gas_u = unyt.unyt_array(readsnap(filename, "u", 0), energy_per_unit_mass)
writer.gas.internal_energy = sw.cosmo_array(
    gas_u.value, gas_u.units, comoving=True, scale_factor=1.0, scale_exponent=-2
)

# We must roll our own smoothing lengths.
n_part = len(gas_coords)
x_range = gas_coords.max() - gas_coords.min()
mean_interparticle_sep = x_range / n_part ** (1 / 3)
writer.gas.smoothing_lengths = sw.cosmo_array(
    np.ones(n_part, dtype=float) * mean_interparticle_sep.value,
    mean_interparticle_sep.units,
    comoving=True,
    scale_factor=1.0,
    scale_exponent=1,
)

writer.gas.particle_ids = sw.cosmo_array(
    readsnap(filename, "pid", 0),
    unyt.dimensionless,
    comoving=True,
    scale_factor=1.0,
    scale_exponent=0,
)


def read_dm_quantity(name, unit, parttype):
    """
    The DM particles are in three sets because of their different masses.
    In SWIFT we have to combine these.
    """
    out = np.concatenate([readsnap(filename, name, p) for p in parttype]) * unit
    return out


for name, parttype in {"dark_matter": [1], "boundary": [2, 3, 5]}.items():
    writer_value = getattr(writer, name)

    dm_coords = read_dm_quantity("pos", length, parttype)
    writer_value.coordinates = sw.cosmo_array(
        dm_coords.value,
        dm_coords.units,
        comoving=True,
        scale_factor=1.0,
        scale_exponent=1,
    )

    dm_vel = read_dm_quantity("vel", velocity, parttype)
    writer_value.velocities = sw.cosmo_array(
        dm_vel.value, dm_vel.units, comoving=True, scale_factor=1.0, scale_exponent=0
    )

    dm_mass = read_dm_quantity("mass", mass, parttype)
    writer_value.masses = sw.cosmo_array(
        dm_mass.value, dm_mass.units, comoving=True, scale_factor=1.0, scale_exponent=0
    )

    dm_pids = read_dm_quantity("pid", unyt.dimensionless, parttype)
    writer_value.particle_ids = sw.cosmo_array(
        dm_pids.value,
        unyt.dimensionless,
        comoving=True,
        scale_factor=1.0,
        scale_exponent=0,
    )

writer.write("nifty.hdf5")
