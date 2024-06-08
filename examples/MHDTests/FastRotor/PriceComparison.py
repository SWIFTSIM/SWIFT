import numpy as np
import h5py
import argparse
import unyt
import matplotlib.pyplot as plt

from swiftsimio import load
from swiftsimio.visualisation.slice import slice_gas

# Parse command line arguments
argparser = argparse.ArgumentParser()
argparser.add_argument("input")
argparser.add_argument("output")
args = argparser.parse_args()

# Load snapshot
filename = args.input
data = load(filename)

# Retrieve particle attributes of interest
rho = data.gas.densities

P = data.gas.pressures

cs = np.sqrt(7.0 * P / (5.0 * rho))

v = data.gas.velocities
normv = np.sqrt(v[:, 0] ** 2 + v[:, 1] ** 2 + v[:, 2] ** 2)

B = data.gas.magnetic_flux_densities

Pmag = 0.5 * (B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2)

# Generate mass weighted maps of quantities of interest
data.gas.mass_weighted_densities = data.gas.masses * rho

data.gas.mass_weighted_pressures = data.gas.masses * P

data.gas.mass_weighted_machnum = data.gas.masses * normv  / cs

data.gas.mass_weighted_Pmag = data.gas.masses * Pmag

common_arguments = dict(data=data, z_slice=0.5 * data.metadata.boxsize[2], resolution=512, parallel=True)

mass_map = slice_gas(**common_arguments, project="masses")

mass_weighted_density_map = slice_gas(
    **common_arguments, project="mass_weighted_densities"
)

mass_weighted_pressure_map = slice_gas(
    **common_arguments, project="mass_weighted_pressures"
)

mass_weighted_machnum_map = slice_gas(**common_arguments, project="mass_weighted_machnum")

mass_weighted_Pmag_map = slice_gas(**common_arguments, project="mass_weighted_Pmag")

# Take out mass dependence
density_map = mass_weighted_density_map / mass_map
pressure_map = mass_weighted_pressure_map / mass_map
machnum_map = mass_weighted_machnum_map / mass_map
Pmag_map = mass_weighted_Pmag_map / mass_map

# Plot maps
plt.rcParams.update({"font.size": 16})
fig, ax = plt.subplots(2, 2, figsize=(12,12))

levels = 30

a00 = ax[0, 0].contour(
    density_map.value.T, levels=levels, linewidths=1.0, colors='k',
)
a01 = ax[0, 1].contour(
    pressure_map.value.T, levels=levels, linewidths=1.0, colors='k',
)
a10 = ax[1, 0].contour(
    machnum_map.value.T, levels=levels, linewidths=1.0, colors='k',
)
a11 = ax[1, 1].contour(
    Pmag_map.value.T, levels=levels, linewidths=1.0, colors='k',
)

for axi in ax:
    for axii in axi:
        axii.set_xticks([])
        axii.set_yticks([])
        axii.set_aspect("equal")

plt.subplots_adjust(wspace=0, hspace=0)

plt.savefig(args.output)
