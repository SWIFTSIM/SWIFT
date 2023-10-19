from swiftsimio import load
from swiftsimio.visualisation.slice import slice_gas
import numpy as np
import sys

filename = sys.argv[1]
data = load(filename)

# First create a mass-weighted temperature dataset
B = data.gas.magnetic_flux_densities
data.gas.mass_weighted_magnetic_pressures = (
    data.gas.masses * (B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2) / 2
)

# Then create a mass-weighted speed dataset
v = data.gas.velocities
data.gas.mass_weighted_speeds = data.gas.masses * (
    v[:, 0] ** 2 + v[:, 1] ** 2 + v[:, 2] ** 2
)

# Map in mass per area
mass_map = slice_gas(
    data,
    z_slice=0.5 * data.metadata.boxsize[2],
    resolution=256,
    project="masses",
    parallel=True,
)

# Map in magnetism squared times mass per area
mass_weighted_magnetic_pressure_map = slice_gas(
    data,
    z_slice=0.5 * data.metadata.boxsize[2],
    resolution=256,
    project="mass_weighted_magnetic_pressures",
    parallel=True,
)
# Map in speed squared times mass per area
mass_weighted_speeds_map = slice_gas(
    data,
    z_slice=0.5 * data.metadata.boxsize[2],
    resolution=256,
    project="mass_weighted_speeds",
    parallel=True,
)

magnetic_pressure_map = mass_weighted_magnetic_pressure_map / mass_map
speed_map = mass_weighted_speeds_map / mass_map

from matplotlib.pyplot import imsave
from matplotlib.colors import LogNorm

# Normalize and save
imsave(sys.argv[2], magnetic_pressure_map.value, cmap="viridis")
imsave(sys.argv[3], speed_map.value, cmap="viridis")
