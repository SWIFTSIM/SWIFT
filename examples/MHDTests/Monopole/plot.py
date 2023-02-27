from swiftsimio import load
from swiftsimio.visualisation.slice import slice_gas
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import ticker

filename = sys.argv[1]
data = load(filename)

print(data.metadata.gas_properties.field_names)

#pos = data.gas.coordinates

# First create a mass-weighted temperature dataset
B        = data.gas.magnetic_flux_density
B_norm_2 = B[:,0]**2 + B[:,1]**2 + B[:,2]**2
P_mag    = B_norm_2 / 2 

h    = data.gas.smoothing_lengths
divB = data.gas.magnetic_divergence 

#error = -10. * np.ones(len(h)) 
#error[B_norm_2 > 0.] = np.log10(np.abs(divB[B_norm_2 > 0.]) * h[B_norm_2 > 0.] / np.sqrt(B_norm_2[B_norm_2 > 0.]))

error = np.log10(divB * h / np.sqrt(B_norm_2))

data.gas.mass_weighted_magnetic_pressures = data.gas.masses * P_mag
data.gas.mass_weighted_error              = data.gas.masses * error

# Then create a mass-weighted pressure dataset
data.gas.mass_weighted_pressures = data.gas.masses * data.gas.pressures

# Then create a mass-weighted speed dataset
v = data.gas.velocities
data.gas.mass_weighted_speeds_x = data.gas.masses * v[:,0]
data.gas.mass_weighted_speeds_y = data.gas.masses * v[:,1]
data.gas.mass_weighted_speeds = data.gas.masses * np.sqrt(v[:,0]**2 + v[:,1]**2 + v[:,2]**2)

# Then create a mass-weighted pressure dataset
data.gas.mass_weighted_densities = data.gas.masses * data.gas.densities

# Map in mass per area
mass_map = slice_gas(data, z_slice=0.5 * data.metadata.boxsize[2], resolution=1024, project="masses", parallel=True)

# Map in density per area
rho_pa_map = slice_gas(data, z_slice=0.5 * data.metadata.boxsize[2], resolution=1024, project="mass_weighted_densities", parallel=True)

# Map in magnetism squared times mass per area
mass_weighted_magnetic_pressure_map = slice_gas(
    data,
    z_slice=0.5 * data.metadata.boxsize[2],
    resolution=1024,
    project="mass_weighted_magnetic_pressures",
    parallel=True
)

# Map in pressure times mass per area
mass_weighted_pressure_map = slice_gas(
    data,
    z_slice=0.5 * data.metadata.boxsize[2],
    resolution=1024,
    project="mass_weighted_pressures",
    parallel=True
)

# Map in speed squared times mass per area
mass_weighted_error_map = slice_gas(
    data,
    z_slice=0.5 * data.metadata.boxsize[2],
    resolution=1024,
    project="mass_weighted_error",
    parallel=True
)

magnetic_pressure_map = mass_weighted_magnetic_pressure_map / mass_map
error_map = mass_weighted_error_map / mass_map
#speed_map = mass_weighted_speeds_map / mass_map
rho_map = rho_pa_map / mass_map
pressure_map = mass_weighted_pressure_map / mass_map

from matplotlib.pyplot import imsave
#from matplotlib.colors import LogNorm

# Normalize and save
#imsave(sys.argv[2], np.rot90(magnetic_pressure_map.value), cmap="jet", vmin=0.0177, vmax=2.642)
#imsave(sys.argv[3], speed_map.value, cmap="viridis")


L = 29

fig, axs = plt.subplots(2, 2, figsize=(12, 10))

"""

im = axs[0,0].contour(np.rot90(rho_map.value), levels, colors='k', vmin=0.483, vmax=12.95)
#plt.colorbar(im, ax = axs[0,0])
im = axs[0,1].contour(np.rot90(pressure_map.value), levels, colors='k', vmin=0.0202, vmax=2.008)
#plt.colorbar(im, ax = axs[0,1])
im = axs[1,0].contour(np.rot90(speed_map.value), levels, colors='k', vmin=0., vmax=2.6) #, vmin=0., vmax=8.18)
#cb = plt.colorbar(im, ax = axs[1,0])
im = axs[1,1].contour(np.rot90(magnetic_pressure_map.value), levels, colors='k', vmin=0.0177, vmax=2.642)
#plt.colorbar(im, ax = axs[1,1])

"""

levels = np.linspace(0.8, 1.2, L)
im = axs[0,0].contourf(np.rot90(rho_map.value), levels=levels, cmap='viridis')
cb = plt.colorbar(im, ax = axs[0,0])

#tick_locator = ticker.MaxNLocator(nbins=10)
#cb.locator = tick_locator
#cb.update_ticks()

levels = np.linspace(0.8, 1.2, L)
im = axs[0,1].contourf(np.rot90(pressure_map.value), levels=levels, cmap='viridis')
cb = plt.colorbar(im, ax = axs[0,1])

#tick_locator = ticker.MaxNLocator(nbins=10)
#cb.locator = tick_locator
#cb.update_ticks()

levels = np.linspace(-1.5, 0.25, L)
im = axs[1,0].contourf(np.rot90(error_map.value), levels=levels, cmap='viridis')
cb = plt.colorbar(im, ax = axs[1,0])

#tick_locator = ticker.MaxNLocator(nbins=10)
#cb.locator = tick_locator
#cb.update_ticks()

levels = np.linspace(0., 60, L)
im = axs[1,1].contourf(np.rot90(magnetic_pressure_map.value), levels=levels, cmap='viridis')
cb = plt.colorbar(im, ax = axs[1,1])

#tick_locator = ticker.MaxNLocator(nbins=10)
#cb.locator = tick_locator
#cb.update_ticks()

#plt.setp(plt.gcf().get_axes(), xticks=[], yticks=[]);

plt.savefig(sys.argv[2])
