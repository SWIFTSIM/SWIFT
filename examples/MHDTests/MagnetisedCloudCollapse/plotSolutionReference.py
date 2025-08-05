# colormaps for the script were chosen to match 1909.09650v2

import numpy as np
import h5py
import sys
import unyt
import matplotlib.pyplot as plt

from swiftsimio import load
from swiftsimio.visualisation.rotation import rotation_matrix_from_vector
from swiftsimio.visualisation.slice import slice_gas

filename = sys.argv[1]
data = load(filename)
center = 0.5 * data.metadata.boxsize

# Constants
G = 6.67430e-11 * unyt.m ** 3 / unyt.kg / unyt.s ** 2  # 6.67430e-8
G = G.to(unyt.cm ** 3 / unyt.g / unyt.s ** 2)
mu0 = 1.25663706127e-1 * unyt.g * unyt.cm / (unyt.s ** 2 * unyt.A ** 2)

# Parameters (taken from Hopkins 2016)
R0 = 4.628516371e16 * unyt.cm
M = 1.99e33 * unyt.g  # total mass of the sphere
rhocloud0 = M / (4 / 3 * np.pi * R0 ** 3)
rhouniform = rhocloud0 / 360
t_ff = np.sqrt(3 / (2 * np.pi * G * rhocloud0)) * unyt.s
tsim = data.metadata.time
toplot = tsim / t_ff
print("Showing results at %3f free fall times" % toplot)

# First create a mass-weighted temperature dataset
r = data.gas.coordinates - center
v = data.gas.velocities
norm_r = np.sqrt(r[:, 0] ** 2 + r[:, 1] ** 2 + r[:, 2] ** 2)
vr = np.sum(v * r, axis=1) / norm_r
B = data.gas.magnetic_flux_densities
J = data.gas.magnetic_flux_curl
normB = np.sqrt(B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2)
normJ = np.sqrt(J[:, 0] ** 2 + J[:, 1] ** 2 + J[:, 2] ** 2)
divB = data.gas.magnetic_divergences
h = data.gas.smoothing_lengths
rho = data.gas.densities


data.gas.mass_weighted_densities = data.gas.masses * data.gas.densities
data.gas.mass_weighted_error = data.gas.masses * np.log10(
    np.maximum(h * abs(divB) / normB, 1e-6)
)
data.gas.mass_weighted_vr = data.gas.masses * vr
data.gas.mass_weighted_Bz = data.gas.masses * abs(B[:, 2])
data.gas.mass_weighted_J = (
    data.gas.masses * np.sqrt(J[:, 0] ** 2 + J[:, 1] ** 2 + J[:, 2] ** 2) / mu0
)

"""
rotation_matrix = [[1,0,0],
                   [0,1,0],
                   [0,0,1]]
"""
rotation_matrix = [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]

Lslice_AU = 1000
Lslice = Lslice_AU * 1.49597871 * 10 ** 13 * unyt.cm

visualise_region_xy = [
    center[0] - Lslice,
    center[0] + Lslice,
    center[1] - Lslice,
    center[1] + Lslice,
]

visualise_region = [
    center[0] - Lslice,
    center[0] + Lslice,
    center[2] - Lslice,
    center[2] + Lslice,
]


common_arguments_xy = dict(
    data=data,
    resolution=512,
    parallel=True,
    region=visualise_region_xy,
    z_slice=center[2],  # 0.0 * unyt.cm,
)

common_arguments = dict(
    data=data,
    resolution=512,
    parallel=True,
    rotation_center=unyt.unyt_array(center),
    rotation_matrix=rotation_matrix,
    region=visualise_region,
    z_slice=0.0 * unyt.cm,
)


# Map in msun / mpc^3

mass_map_xy = slice_gas(**common_arguments_xy, project="masses")
mass_weighted_density_map_xy = slice_gas(
    **common_arguments_xy, project="mass_weighted_densities"
)

mass_map = slice_gas(**common_arguments, project="masses")
mass_weighted_density_map = slice_gas(
    **common_arguments, project="mass_weighted_densities"
)
mass_weighted_error_map = slice_gas(**common_arguments, project="mass_weighted_error")
mass_weighted_vr_map = slice_gas(**common_arguments, project="mass_weighted_vr")
mass_weighted_Bz_map = slice_gas(**common_arguments, project="mass_weighted_Bz")
mass_weighted_J_map = slice_gas(**common_arguments, project="mass_weighted_J")

density_map_xy = mass_weighted_density_map_xy / mass_map_xy
density_map = mass_weighted_density_map / mass_map
error_map = mass_weighted_error_map / mass_map
vr_map = mass_weighted_vr_map / mass_map
Bz_map = mass_weighted_Bz_map / mass_map
J_map = mass_weighted_J_map / mass_map

density_map_xy.convert_to_units(unyt.g * unyt.cm ** (-3))
density_map.convert_to_units(unyt.g * unyt.cm ** (-3))
vr_map.convert_to_units(unyt.km / unyt.s)
Bz_map.convert_to_units(1e-7 * unyt.g / (unyt.A * unyt.s * unyt.s))
J_map.convert_to_units(unyt.A / (unyt.m ** 2))

fig, ax = plt.subplots(3, 2, sharey=True, figsize=(10, 15))

fig.suptitle("Plotting at %.3f free fall times" % toplot)

a00 = ax[0, 0].contourf(
    np.log10(density_map.value),
    cmap="RdBu",
    extend="both",
    levels=np.linspace(-17.0, -11.0, 100),
)

a01 = ax[0, 1].contourf(
    np.log10(density_map_xy.value),
    cmap="RdBu",
    extend="both",
    levels=np.linspace(-17.0, -11.0, 100),
)

a10 = ax[1, 0].contourf(
    vr_map.value, cmap="CMRmap", extend="both", levels=np.linspace(-1.0, 1.0, 100)
)
a11 = ax[1, 1].contourf(
    error_map.value, cmap="jet", extend="both", levels=np.arange(-6.0, 3.0, 1.0)
)
a20 = ax[2, 0].contourf(
    np.log10(np.maximum(J_map.value, -15)),
    cmap="nipy_spectral",
    extend="both",
    levels=np.linspace(-15, -11, 100),
)
a21 = ax[2, 1].contourf(
    np.log10(np.maximum(Bz_map.value, -6.9)),
    cmap="nipy_spectral_r",
    extend="both",
    levels=np.linspace(2.0, 5.0, 100),
)

locs = [512 / 4, 512 / 2, 3 * 512 / 4]
labels = [-Lslice_AU / 2, 0, Lslice_AU / 2]

for ii in range(3):
    ax[ii, 0].set_ylabel(r"$z$ [A.U.]")
    ax[ii, 0].set_yticks(locs, labels)

for ii in range(3):
    for jj in range(2):
        ax[ii, jj].set_xlabel(r"$x$ [A.U.]")
        ax[ii, jj].set_xticks(locs, labels)
        ax[ii, jj].set_aspect("equal")

cbar1 = plt.colorbar(
    a00, ax=ax[0, 0], fraction=0.046, pad=0.04, ticks=np.linspace(-17, -11, 7)
)
cbar1.set_label(r"$\mathrm{log}_{10} \rho \; [g {cm}^{-3}]$")
cbar2 = plt.colorbar(
    a01, ax=ax[0, 1], fraction=0.046, pad=0.04, ticks=np.linspace(-17, -11, 7)
)
ax[0, 1].set_ylabel(r"$y$ [A.U.]")
cbar2.set_label(r"$\mathrm{log}_{10} \rho \; [g {cm}^{-3}]$")
cbar3 = plt.colorbar(
    a10, ax=ax[1, 0], fraction=0.046, pad=0.04, ticks=np.linspace(-1.0, 1.0, 9)
)
cbar3.set_label(r"$\mathrm{log}_{10} v_r \; [km s^{-1}]$")

cbar4 = plt.colorbar(a11, ax=ax[1, 1], fraction=0.046, pad=0.04)
cbar4.set_label(r"$\mathrm{log}_{10} \frac{h \: \nabla \cdot B}{|B|}$")
cbar5 = plt.colorbar(
    a20, ax=ax[2, 0], fraction=0.046, pad=0.04, ticks=np.linspace(-15, -11, 5)
)
cbar5.set_label(r"$\mathrm{log}_{10} |J| \; [A/m^2]$")
cbar6 = plt.colorbar(
    a21, ax=ax[2, 1], fraction=0.046, pad=0.04, ticks=np.linspace(2, 5, 4)
)
cbar6.set_label(r"$\mathrm{log}_{10} B_z \; [\mu G]$")

plt.tight_layout()

plt.savefig(sys.argv[2])
