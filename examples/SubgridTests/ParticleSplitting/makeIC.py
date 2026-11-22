"""
Creates an IC for the particle splitting test.

This test creates a single heavy particle, that should
split into 16 counterparts, surrounded by a grid of
16 * 16 * 16 particles.
"""

import numpy as np
import unyt

import swiftsimio as sw

NUMBER_OF_PARTICLES = 16  # in 1D
TOTAL_NUMBER_OF_PARTICLES = NUMBER_OF_PARTICLES**3
PARTICLE_OVER_MASS_FACTOR = 16  # how many times over-massive is the main particle?
TIME_STEP = 0.1  # in seconds

boxsize = sw.cosmo_array(
    [1, 1, 1], unyt.cm, comoving=True, scale_factor=1.0, scale_exponent=1
)
dataset = sw.Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    boxsize=boxsize,
    dimension=3,
)

base_coordinates = np.arange(0.0, 1.0, 1.0 / NUMBER_OF_PARTICLES)
dataset.gas.coordinates = sw.cosmo_array(
    np.array([x.flatten() for x in np.meshgrid(*[base_coordinates] * 3)]).T,
    unyt.cm,
    comoving=True,
    scale_factor=1.0,
    scale_exponent=1,
)
dataset.gas.generate_smoothing_lengths()

base_velocities = np.zeros_like(base_coordinates)
dataset.gas.velocities = sw.cosmo_array(
    np.array([x.flatten() for x in np.meshgrid(*[base_velocities] * 3)]).T,
    unyt.cm / unyt.s,
    comoving=True,
    scale_factor=1.0,
    scale_exponent=0,
)

# Set the particle with the highest mass to be in the centre of
# the volume for easier plotting later
special_particle = (
    np.linalg.norm(
        dataset.gas.coordinates
        - sw.cosmo_array(
            [0.5, 0.5, 0.5], unyt.cm, comoving=True, scale_factor=1.0, scale_exponent=1
        ),
        axis=1,
    )
).argmin()

base_masses = np.ones(TOTAL_NUMBER_OF_PARTICLES, dtype=np.float32)
base_masses[special_particle] = np.float32(PARTICLE_OVER_MASS_FACTOR)
dataset.gas.masses = sw.cosmo_array(
    base_masses, unyt.g, comoving=True, scale_factor=1.0, scale_exponent=0
)

# Set internal energy to be consistent with a CFL time-step of the
# required length
h0 = dataset.gas.smoothing_lengths[0].to_value(unyt.cm)
internal_energy = (h0 * 0.1 / TIME_STEP) ** 2 / (5 / 3 * (5 / 3 - 1))
internal_energies = (
    np.ones(TOTAL_NUMBER_OF_PARTICLES, dtype=np.float32) * internal_energy
)
dataset.gas.internal_energy = sw.cosmo_array(
    internal_energies,
    unyt.cm**2 / unyt.s**2,
    comoving=True,
    scale_factor=1.0,
    scale_exponent=-2,
)

dataset.gas.particle_ids = sw.cosmo_array(
    1 + np.arange(TOTAL_NUMBER_OF_PARTICLES),
    unyt.dimensionless,
    comoving=True,
    scale_factor=1.0,
    scale_exponent=0,
)

dataset.write("particleSplitting.hdf5")
