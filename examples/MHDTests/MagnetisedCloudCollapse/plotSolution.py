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

R0 = 0.015 * 3.086e18 * unyt.cm

tff = 3e4 * 3.156e7 * unyt.s
tsim = data.metadata.time
tplot = tsim / tff
print("Showing results at %2f free fall times" % tplot)

# First create a mass-weighted temperature dataset
r = data.gas.coordinates - center
v = data.gas.velocities
norm_r = np.sqrt(r[:, 0] ** 2 + r[:, 1] ** 2 + r[:, 2] ** 2)
vr = np.sum(v * r, axis=1) / norm_r

B = data.gas.magnetic_flux_density
normB = np.sqrt(B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2)
divB = data.gas.magnetic_divergence
h = data.gas.smoothing_lengths

data.gas.mass_weighted_densities = data.gas.masses * data.gas.densities
data.gas.mass_weighted_error = data.gas.masses * np.log10(
    np.maximum(h * abs(divB) / normB, 1e-6)
)
data.gas.mass_weighted_vr = (
    data.gas.masses * vr
)  # np.sqrt(vr[:,0]**2 + vr[:,1]**2 + vr[:,2]**2)
data.gas.mass_weighted_Btheta = data.gas.masses * np.sqrt(B[:, 0] ** 2 + B[:, 1] ** 2)
data.gas.mass_weighted_Bz = data.gas.masses * abs(B[:, 2])
"""
rotation_matrix = [[1,0,0],
                   [0,1,0],
                   [0,0,1]]
"""
rotation_matrix = [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]

visualise_region = [
    center[0] - 0.15 * R0,
    center[0] + 0.15 * R0,
    center[2] - 0.15 * R0,
    center[2] + 0.15 * R0,
]

visualise_region_zoom = [
    center[0] - 0.015 * R0,
    center[0] + 0.015 * R0,
    center[2] - 0.015 * R0,
    center[2] + 0.015 * R0,
]

"""
visualise_region = [
   center[0] - 2.0*R0, center[0] + 2.0*R0,
   center[2] - 2.0*R0, center[2] + 2.0*R0
]

visualise_region_zoom = [
   center[0] - 0.2*R0, center[0] + 0.2*R0,
   center[2] - 0.2*R0, center[2] + 0.2*R0
]
"""

common_arguments = dict(
    data=data,
    resolution=512,
    parallel=True,
    rotation_center=unyt.unyt_array(center),
    rotation_matrix=rotation_matrix,
    region=visualise_region,
    z_slice=0.01 * unyt.cm,
)
common_arguments_zoom = dict(
    data=data,
    resolution=512,
    parallel=True,
    rotation_center=unyt.unyt_array(center),
    rotation_matrix=rotation_matrix,
    region=visualise_region_zoom,
    z_slice=0.01 * unyt.cm,
)

# Map in msun / mpc^3
mass_map = slice_gas(**common_arguments, project="masses")

mass_zoom_map = slice_gas(**common_arguments_zoom, project="masses")

mass_weighted_density_map = slice_gas(
    **common_arguments, project="mass_weighted_densities"
)

mass_weighted_error_map = slice_gas(**common_arguments, project="mass_weighted_error")

mass_weighted_vr_map = slice_gas(**common_arguments, project="mass_weighted_vr")

mass_weighted_vr_zoom_map = slice_gas(
    **common_arguments_zoom, project="mass_weighted_vr"
)

mass_weighted_Btheta_map = slice_gas(**common_arguments, project="mass_weighted_Btheta")

mass_weighted_Bz_map = slice_gas(**common_arguments, project="mass_weighted_Bz")

density_map = mass_weighted_density_map / mass_map
error_map = mass_weighted_error_map / mass_map
vr_map = mass_weighted_vr_map / mass_map

# vr_map[vr_map.value == -np.inf] = 0.0 * unyt.cm / unyt.s

vr_zoom_map = mass_weighted_vr_zoom_map / mass_zoom_map
Btheta_map = mass_weighted_Btheta_map / mass_map
Bz_map = mass_weighted_Bz_map / mass_map

density_map.convert_to_units(unyt.g * unyt.cm ** (-3))
vr_map.convert_to_units(unyt.km / unyt.s)
vr_zoom_map.convert_to_units(unyt.km / unyt.s)
Btheta_map.convert_to_units(1e-7 * unyt.g / (unyt.statA * unyt.s * unyt.s))
Bz_map.convert_to_units(1e-7 * unyt.g / (unyt.statA * unyt.s * unyt.s))

fig, ax = plt.subplots(3, 2, sharey=True, figsize=(10, 15))

fig.suptitle("Plotting at %.2f free fall times" % tplot)

a00 = ax[0, 0].contourf(
    np.log10(density_map.value),
    cmap="jet",
    extend="both",
    levels=np.linspace(-17.0, -11.0, 100),
)
a01 = ax[0, 1].contourf(
    error_map.value, cmap="jet", extend="both", levels=np.arange(-6.0, 3.0, 1.0)
)
a10 = ax[1, 0].contourf(
    vr_map.value, cmap="jet", extend="both", levels=np.linspace(-1.5, 2.5, 100)
)
a11 = ax[1, 1].contourf(
    vr_zoom_map.value, cmap="jet", extend="both", levels=np.linspace(-1.5, 2.5, 100)
)
a20 = ax[2, 0].contourf(
    np.log10(np.maximum(Btheta_map.value, 1.0)),
    cmap="jet",
    extend="both",
    levels=np.linspace(0.0, 5.0, 100),
)
a21 = ax[2, 1].contourf(
    np.log10(np.maximum(Bz_map.value, 1.0)),
    cmap="jet",
    extend="both",
    levels=np.linspace(1.0, 4.5, 100),
)

# locs = [10.67, 32, 53.33]
# locs = [21.33, 64, 106.67]
locs = [85.33, 256.00, 426.67]
# locs   = [170.67, 512.00, 853.33]
# locs   = [341.33, 1024, 1706.67]
labels = [-0.1, 0.0, 0.1]

# locs   = [128, 256, 384]
# labels = [-1.0, 0.0, 1.0]

for ii in range(3):
    ax[ii, 0].set_ylabel(r"$z/R_0$")
    ax[ii, 0].set_yticks(locs, labels)

for ii in range(3):
    for jj in range(2):
        # R = 0.001/15 * 1707
        # print(R)
        # circle = plt.Circle((256,256),R, color='k', fill=False)
        # ax[ii,jj].add_patch(circle)
        ax[ii, jj].set_xlabel(r"$x/R_0$")
        ax[ii, jj].set_xticks(locs, labels)
        ax[ii, jj].set_aspect("equal")

cbar1 = plt.colorbar(
    a00, ax=ax[0, 0], fraction=0.046, pad=0.04, ticks=np.linspace(-17, -11, 7)
)
cbar1.set_label(r"$\mathrm{log}_{10} \rho \; [g {cm}^{-3}]$")
cbar2 = plt.colorbar(a01, ax=ax[0, 1], fraction=0.046, pad=0.04)
cbar2.set_label(r"$\mathrm{log}_{10} \frac{h \: \nabla \cdot B}{|B|}$")
cbar3 = plt.colorbar(
    a10, ax=ax[1, 0], fraction=0.046, pad=0.04, ticks=np.linspace(-1.5, 2.5, 9)
)
cbar3.set_label(r"$v_r \; [km s^{-1}$")
cbar4 = plt.colorbar(
    a11, ax=ax[1, 1], fraction=0.046, pad=0.04, ticks=np.linspace(-1.5, 2.5, 9)
)
cbar4.set_label(r"$v_r \; [km s^{-1}$")
cbar5 = plt.colorbar(
    a20, ax=ax[2, 0], fraction=0.046, pad=0.04, ticks=np.linspace(0.0, 5.0, 6)
)
cbar5.set_label(r"$\mathrm{log}_{10} B_\theta \; [\mu G]$")
cbar6 = plt.colorbar(
    a21, ax=ax[2, 1], fraction=0.046, pad=0.04, ticks=np.linspace(1.0, 4.5, 8)
)
cbar6.set_label(r"$\mathrm{log}_{10} B_z \; [\mu G]$")

plt.tight_layout()

plt.savefig(sys.argv[2])
