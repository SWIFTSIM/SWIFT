from swiftsimio import load
from swiftsimio.visualisation.projection import project_gas
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
# v = data.gas.velocities
# data.gas.mass_weighted_speeds = data.gas.masses * (v[:,0]**2 + v[:,1]**2 + v[:,2]**2)

# Map in mass per area
mass_map = project_gas(data, resolution=1024, project="masses", parallel=True)
# Map in magnetism squared times mass per area
mass_weighted_magnetic_pressure_map = project_gas(
    data, resolution=1024, project="mass_weighted_magnetic_pressures", parallel=True
)
"""
# Map in speed squared times mass per area
mass_weighted_speeds_map = project_gas(
    data,
    resolution=1024,
    project="mass_weighted_speeds",
    parallel=True
)
"""
magnetic_pressure_map = mass_weighted_magnetic_pressure_map / mass_map
# speed_map = mass_weighted_speeds_map / mass_map

from matplotlib.pyplot import imsave

# from matplotlib.colors import LogNorm

# Normalize and save
imsave(sys.argv[2], magnetic_pressure_map.value.T, cmap="viridis")
# imsave(sys.argv[3], speed_map.value, cmap="viridis")
