from swiftsimio import load
from swiftsimio.visualisation.projection import project_gas
import numpy as np
import sys

input_filename_base = "MagneticBlastWave_LR_"

for ii in range(41):
	
	print(ii)

	filename = input_filename_base + str(ii).zfill(4) + ".hdf5"
	data = load(filename)

	# First create a mass-weighted temperature dataset
	rho = data.gas.densities
	data.gas.mass_weighted_densities = data.gas.masses * rho

	# Map in mass per area
	mass_map = project_gas(data, resolution=1024, project="masses", parallel=True)
	# Map in magnetism squared times mass per area
	mass_weighted_density_map = project_gas(
    	data,
    	resolution=1024,
    	project="mass_weighted_densities",
    	parallel=True
	)

	density_map = mass_weighted_density_map / mass_map

	from matplotlib.pyplot import imsave

	# Normalize and save
	imsave(input_filename_base + str(ii).zfill(4) + ".png", np.rot90(density_map.value), cmap="jet", vmin=0.1, vmax = 2.4)
