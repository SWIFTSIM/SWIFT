"""
Creates an IC for the particle splitting test.

This test creates a single heavy particle, that should
split into 16 counterparts, surrounded by a grid of
16 * 16 * 16 particles.
"""

import numpy as np
import unyt

from swiftsimio import Writer

NUMBER_OF_PARTICLES = 16  # in 1D
TOTAL_NUMBER_OF_PARTICLES = NUMBER_OF_PARTICLES ** 3
PARTICLE_OVER_MASS_FACTOR = 16  # how many times over-massive is the main particle?
TIME_STEP = 0.1  # in seconds

dataset = Writer(
    unit_system=unyt.unit_systems.cgs_unit_system,
    box_size=unyt.unyt_array([1, 1, 1], "cm"),
    dimension=3,
)

base_coordinates = np.arange(0.0, 1.0, 1.0 / NUMBER_OF_PARTICLES)
dataset.gas.coordinates = unyt.unyt_array(
    [x.flatten() for x in np.meshgrid(*[base_coordinates] * 3)], "cm"
).T
dataset.gas.generate_smoothing_lengths(boxsize=dataset.box_size, dimension=3)

base_velocities = np.zeros_like(base_coordinates)
dataset.gas.velocities = unyt.unyt_array(
    [x.flatten() for x in np.meshgrid(*[base_velocities] * 3)], "cm/s"
).T

# Set the particle with the highest mass to be in the centre of
# the volume for easier plotting later
special_particle = (
    np.linalg.norm(dataset.gas.coordinates - unyt.unyt_quantity(0.5, "cm"), axis=1)
).argmin()

base_masses = np.ones(TOTAL_NUMBER_OF_PARTICLES, dtype=np.float32)
base_masses[special_particle] = np.float32(PARTICLE_OVER_MASS_FACTOR)
dataset.gas.masses = unyt.unyt_array(base_masses, "g")

# Set internal energy to be consistent with a CFL time-step of the
# required length
internal_energy = (dataset.gas.smoothing_length[0] * 0.1 / TIME_STEP) ** 2 / (
    5 / 3 * (5 / 3 - 1)
)
internal_energies = (
    np.ones(TOTAL_NUMBER_OF_PARTICLES, dtype=np.float32) * internal_energy
)
dataset.gas.internal_energy = unyt.unyt_array(internal_energies, "cm**2 / s**2")

dataset.gas.particle_ids = unyt.unyt_array(
    np.arange(TOTAL_NUMBER_OF_PARTICLES), "dimensionless"
)

dataset.write("particleSplitting.hdf5")

