from swiftsimio import load
from swiftsimio.visualisation.projection import project_gas
import numpy as np
import sys
import matplotlib.pyplot as plt
from matplotlib import ticker

filename = sys.argv[1]
data = load(filename)

print(data.metadata.gas_properties.field_names)

divB                        = data.gas.magnetic_divergence 
data.gas.mass_weighted_divB = data.gas.masses * divB

# Map in mass per area
mass_map = project_gas(data, resolution=1024, project="masses", parallel=True)

# Map in speed squared times mass per area
mass_weighted_divB_map = project_gas(
    data,
    resolution=1024,
    project="mass_weighted_divB",
    parallel=True
)

divB_map = mass_weighted_divB_map / mass_map

from matplotlib.pyplot import imsave

imsave(sys.argv[2], np.rot90(divB_map.value), cmap="gray", vmin=-1.5, vmax=1.5)

quit()

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
