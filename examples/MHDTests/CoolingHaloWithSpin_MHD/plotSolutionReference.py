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

time = data.metadata.time
time.convert_to_units(unyt.Gyr)
print("Plotting at %2f Gyr" % time)

# Retrieve some information about the simulation run
artDiffusion = data.metadata.hydro_scheme["Artificial Diffusion Constant"]
dedHyp = data.metadata.hydro_scheme["Dedner Hyperbolic Constant"]
dedHypDivv = data.metadata.hydro_scheme["Dedner Hyperbolic div(v) Constant"]
dedPar = data.metadata.hydro_scheme["Dedner Parabolic Constant"]
eta = data.metadata.hydro_scheme["Resistive Eta"]
git = data.metadata.code["Git Revision"]
gitBranch = data.metadata.code["Git Branch"]
hydroScheme = data.metadata.hydro_scheme["Scheme"]
kernel = data.metadata.hydro_scheme["Kernel function"]
neighbours = data.metadata.hydro_scheme["Kernel target N_ngb"]


# Constants
G = 6.67430e-11 * unyt.m ** 3 / unyt.kg / unyt.s ** 2  # 6.67430e-8
G = G.to(unyt.cm ** 3 / unyt.g / unyt.s ** 2)
mu0 = 1.25663706127e-1 * unyt.g * unyt.cm / (unyt.s ** 2 * unyt.statA ** 2)

# First create a mass-weighted temperature dataset
r = data.gas.coordinates - center
v = data.gas.velocities
norm_r = np.sqrt(r[:, 0] ** 2 + r[:, 1] ** 2 + r[:, 2] ** 2)
vr = np.sum(v * r, axis=1) / norm_r
normv = np.sqrt(v[:, 0] ** 2 + v[:, 1] ** 2 + v[:, 2] ** 2)
B = data.gas.magnetic_flux_densities
# B = v.value * 0.0 * 1e-7*unyt.g / (unyt.statA * unyt.s * unyt.s)
# J = data.gas.magnetic_flux_curl
normB = np.sqrt(B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2)
# normJ = np.sqrt(J[:, 0] ** 2 + J[:, 1] ** 2 + J[:, 2] ** 2)
divB = data.gas.magnetic_divergences
# divB = norm_r.value * 0.0 * 1e-7*unyt.g / (unyt.statA * unyt.s * unyt.s* unyt.cm)
h = data.gas.smoothing_lengths
rho = data.gas.densities
Np = len(h)
print("Number of particles %E" % len(data.gas.masses))
print("Total mass %E Msol" % np.sum(data.gas.masses.to(unyt.Msun).value))
print("Gas particle mass %E Msol" % np.mean(data.gas.masses.to(unyt.Msun).value))

plasmaBeta = data.gas.pressures / (normB ** 2 / (2 * mu0))

data.gas.mass_weighted_densities = data.gas.masses * data.gas.densities
data.gas.mass_weighted_error = data.gas.masses * np.log10(
    np.maximum(h * abs(divB) / normB, 1e-6)
)
data.gas.mass_weighted_vr = data.gas.masses * vr
data.gas.mass_weighted_Bz = data.gas.masses * abs(B[:, 2])
data.gas.mass_weighted_Bx = data.gas.masses * abs(B[:, 0])
data.gas.mass_weighted_By = data.gas.masses * abs(B[:, 1])
# data.gas.mass_weighted_J = (
#    data.gas.masses * np.sqrt(J[:, 0] ** 2 + J[:, 1] ** 2 + J[:, 2] ** 2) / mu0
# )

data.gas.mass_weighted_plasmaBeta = data.gas.masses * plasmaBeta
data.gas.mass_weighted_normv = data.gas.masses * normv
data.gas.mass_weighted_normB = data.gas.masses * normB

"""
rotation_matrix = [[1,0,0],
                   [0,1,0],
                   [0,0,1]]
"""
rotation_matrix = [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]

# set plot area
Lslice_kPc = 20.0  # 2*21.5
Lslice = Lslice_kPc * 3.08e18 * 1e3 * unyt.cm

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
mass_weighted_error_map_xy = slice_gas(
    **common_arguments_xy, project="mass_weighted_error"
)
mass_weighted_vr_map = slice_gas(**common_arguments, project="mass_weighted_vr")
mass_weighted_Bz_map = slice_gas(**common_arguments, project="mass_weighted_Bz")
mass_weighted_Bx_map = slice_gas(**common_arguments, project="mass_weighted_Bx")
mass_weighted_By_map = slice_gas(**common_arguments, project="mass_weighted_By")
# mass_weighted_J_map = slice_gas(**common_arguments, project="mass_weighted_J")

mass_weighted_plasmaBeta_map = slice_gas(
    **common_arguments, project="mass_weighted_plasmaBeta"
)
mass_weighted_normv_map = slice_gas(**common_arguments, project="mass_weighted_normv")
mass_weighted_normv_map_xy = slice_gas(
    **common_arguments_xy, project="mass_weighted_normv"
)
mass_weighted_normB_map = slice_gas(**common_arguments, project="mass_weighted_normB")
mass_weighted_normB_map_xy = slice_gas(
    **common_arguments_xy, project="mass_weighted_normB"
)
mass_weighted_Bz_map_xy = slice_gas(**common_arguments_xy, project="mass_weighted_Bz")
mass_weighted_Bx_map_xy = slice_gas(**common_arguments_xy, project="mass_weighted_Bx")
mass_weighted_By_map_xy = slice_gas(**common_arguments_xy, project="mass_weighted_By")

density_map_xy = mass_weighted_density_map_xy / mass_map_xy
density_map = mass_weighted_density_map / mass_map
error_map = mass_weighted_error_map / mass_map
error_map_xy = mass_weighted_error_map_xy / mass_map_xy
vr_map = mass_weighted_vr_map / mass_map
Bz_map = mass_weighted_Bz_map / mass_map
Bx_map = mass_weighted_Bx_map / mass_map
By_map = mass_weighted_By_map / mass_map
Bz_map_xy = mass_weighted_Bz_map_xy / mass_map_xy
Bx_map_xy = mass_weighted_Bx_map_xy / mass_map_xy
By_map_xy = mass_weighted_By_map_xy / mass_map_xy
# J_map = mass_weighted_J_map / mass_map
plasmaBeta_map = mass_weighted_plasmaBeta_map / mass_map
normv_map = mass_weighted_normv_map / mass_map
normv_map_xy = mass_weighted_normv_map_xy / mass_map_xy
normB_map = mass_weighted_normB_map / mass_map
normB_map_xy = mass_weighted_normB_map_xy / mass_map_xy

# density_map_xy.convert_to_units(unyt.g * unyt.cm ** (-3))
# density_map.convert_to_units(unyt.g * unyt.cm ** (-3))

density_map_xy.convert_to_units(1e10 * (1e33 * unyt.g) * (3.08e21 * unyt.cm) ** (-3))
density_map.convert_to_units(1e10 * (1e33 * unyt.g) * (3.08e21 * unyt.cm) ** (-3))

nH_map_xy = density_map_xy / (1.67e-24 * unyt.g)
nH_map = density_map / (1.67e-24 * unyt.g)

vr_map.convert_to_units(unyt.km / unyt.s)
Bz_map.convert_to_units(1e-7 * unyt.g / (unyt.statA * unyt.s * unyt.s))
Bx_map.convert_to_units(1e-7 * unyt.g / (unyt.statA * unyt.s * unyt.s))
By_map.convert_to_units(1e-7 * unyt.g / (unyt.statA * unyt.s * unyt.s))
Bz_map_xy.convert_to_units(1e-7 * unyt.g / (unyt.statA * unyt.s * unyt.s))
Bx_map_xy.convert_to_units(1e-7 * unyt.g / (unyt.statA * unyt.s * unyt.s))
By_map_xy.convert_to_units(1e-7 * unyt.g / (unyt.statA * unyt.s * unyt.s))
# J_map.convert_to_units(unyt.statA / (unyt.m ** 2))
normv_map.convert_to_units(unyt.km / unyt.s)
normv_map_xy.convert_to_units(unyt.km / unyt.s)
normB_map.convert_to_units(1e-7 * unyt.g / (unyt.statA * unyt.s * unyt.s))
normB_map_xy.convert_to_units(1e-7 * unyt.g / (unyt.statA * unyt.s * unyt.s))
# normB_map = np.sqrt(Bx_map**2+By_map**2+Bz_map**2)
# normB_map_xy = np.sqrt(Bx_map_xy**2+By_map_xy**2+Bz_map_xy**2)

nx = 2
ny = 4
fig, ax = plt.subplots(ny, nx, sharey=True, figsize=(5 * nx, 5 * ny))

# fig.suptitle("Plotting at %.3f free fall times" % toplot)

a00 = ax[0, 0].contourf(
    np.log10(nH_map.value),
    # np.log10(density_map.value),
    cmap="jet",  # "plasma",#"gist_stern",#"RdBu",
    extend="both",
    levels=np.linspace(-4, 2, 100),
    # levels=np.linspace(-7, -1, 100),
)

a01 = ax[0, 1].contourf(
    np.log10(nH_map_xy.value).T,
    # np.log10(density_map_xy.value).T,
    cmap="jet",  # "plasma",#"gist_stern", #"RdBu",
    extend="both",
    levels=np.linspace(-4, 2, 100),
    # levels=np.linspace(-7, -1, 100),
)

# a10 = ax[1, 0].contourf(
#    np.log10(np.maximum(normv_map.value,-1)), cmap="cividis", extend="both", levels=np.linspace(1.0, 2.0, 100)
# )
# a11 = ax[1, 1].contourf(
#    np.log10(np.maximum(normv_map_xy.value,-1)), cmap="cividis", extend="both", levels=np.linspace(1.0, 2.0, 100)
# )
# a11 = ax[1, 1].contourf(
#    error_map.value, cmap="jet", extend="both", levels=np.arange(-6.0, 3.0, 1.0)
# )
# a20 = ax[2, 0].contourf(
#    np.log10(np.maximum(plasmaBeta_map.value, -4)),
#    cmap="gist_stern",
#    extend="both",
#    levels=np.linspace(-3, 4, 100),
# )
a10 = ax[1, 0].contourf(
    np.maximum(np.log10(normB_map.value), -6),
    cmap="jet",  # "viridis", #"nipy_spectral_r",
    extend="both",
    levels=np.linspace(-2, 2, 100),
)
a11 = ax[1, 1].contourf(
    np.maximum(np.log10(normB_map_xy.value), -6).T,
    cmap="jet",  # "viridis", #"nipy_spectral_r",
    extend="both",
    levels=np.linspace(-2, 2, 100),
)

a20 = ax[2, 0].contourf(
    error_map.value, cmap="jet", extend="both", levels=np.arange(-6.0, 3.0, 1.0)
)
a21 = ax[2, 1].contourf(
    (error_map_xy.value).T, cmap="jet", extend="both", levels=np.arange(-6.0, 3.0, 1.0)
)

a30 = ax[3, 0].contourf(
    np.maximum(np.log10(plasmaBeta_map.value), -10),
    cmap="bwr",
    extend="both",
    levels=np.linspace(0.0, 1.0, 100),
)
locs = [512 / 4, 512 / 2, 3 * 512 / 4]
locs = [512 / 4, 512 / 2, 3 * 512 / 4]
labels = [-Lslice_kPc / 2, 0, Lslice_kPc / 2]

for ii in range(ny):
    ax[ii, 0].set_ylabel(r"$z$ [kPc]")
    ax[ii, 0].set_yticks(locs, labels)

for ii in range(ny):
    for jj in range(nx):
        ax[ii, jj].set_xlabel(r"$x$ [kPc]")
        ax[ii, jj].set_xticks(locs, labels)
        # ax[ii, jj].set_aspect("equal")
        ax[ii, jj].set_aspect("equal", adjustable="box")
        ax[ii, jj].set_xlim(0, 511)
        ax[ii, jj].set_ylim(0, 511)


ticks_rho = [-4, -3, -2, -1, 0, 1, 2]
# ticks_rho = [-7,-6,-5,-4,-3,-2,-1]
cbar1 = plt.colorbar(a00, ax=ax[0, 0], fraction=0.046, pad=0.04, ticks=ticks_rho)
cbar1.set_label(r"$\mathrm{log}_{10} \rho / m_H \; [ {cm}^{-3}]$")
# cbar1.set_label(r"$\mathrm{log}_{10} \rho \; [10^{10} {M}_{sol} {kpc}^{-3}]$")
cbar2 = plt.colorbar(a01, ax=ax[0, 1], fraction=0.046, pad=0.04, ticks=ticks_rho)
ax[0, 1].set_ylabel(r"$y$ [kPc]")
cbar2.set_label(r"$\mathrm{log}_{10} \rho / m_H \; [ {cm}^{-3}]$")
# cbar2.set_label(r"$\mathrm{log}_{10} \rho \; [10^{10} {M}_{sol} {kpc}^{-3}]$")

# cbar3 = plt.colorbar(
#    a10, ax=ax[1, 0], fraction=0.046, pad=0.04, ticks=np.linspace(1, 2, 2)
# )
# cbar3.set_label(r"$\mathrm{log}_{10} |v| \; [km s^{-1}]$")
#
# cbar4 = plt.colorbar(
#    a11, ax=ax[1, 1], fraction=0.046, pad=0.04, ticks=np.linspace(1, 2, 2)
# )
# cbar4.set_label(r"$\mathrm{log}_{10} |v| \; [km s^{-1}]$")
ax[1, 1].set_ylabel(r"$y$ [kPc]")
# cbar4 = plt.colorbar(a11, ax=ax[1, 1], fraction=0.046, pad=0.04)
# cbar4.set_label(r"$\mathrm{log}_{10} \frac{h \: \nabla \cdot B}{|B|}$")
# ticks_Beta = [-3,-1.83,-0.667,0.500,1.67,2.83,4.0]
# cbar5 = plt.colorbar(
#    a20, ax=ax[2, 0], fraction=0.046, pad=0.04, ticks=ticks_Beta
# )
# cbar5.set_label(r"$\mathrm{log}_{10} \beta \;$")
ticks_B = [-2, -1, 0, 1, 2]
cbar3 = plt.colorbar(a11, ax=ax[1, 0], fraction=0.046, pad=0.04, ticks=ticks_B)
cbar3.set_label(r"$\mathrm{log}_{10} B \; [\mu G]$")
cbar4 = plt.colorbar(a11, ax=ax[1, 1], fraction=0.046, pad=0.04, ticks=ticks_B)
cbar4.set_label(r"$\mathrm{log}_{10} B \; [\mu G]$")

ax[2, 1].set_ylabel(r"$y$ [kPc]")

cbar5 = plt.colorbar(a20, ax=ax[2, 0], fraction=0.046, pad=0.04)
cbar5.set_label(r"$\mathrm{log}_{10} \frac{h \: \nabla \cdot B}{|B|}$")

cbar6 = plt.colorbar(a21, ax=ax[2, 1], fraction=0.046, pad=0.04)
cbar6.set_label(r"$\mathrm{log}_{10} \frac{h \: \nabla \cdot B}{|B|}$")

ticksBeta = [0, 1]
cbar7 = plt.colorbar(a30, ax=ax[3, 0], fraction=0.046, pad=0.04, ticks=ticksBeta)
cbar7.set_label(r"$\mathrm{log}_{10} \beta $")

# add streamlines

data.gas.mass_weighted_vx = data.gas.masses * v[:, 0]
data.gas.mass_weighted_vy = data.gas.masses * v[:, 1]
mass_weighted_vx_map_xy = slice_gas(**common_arguments_xy, project="mass_weighted_vx")
mass_weighted_vy_map_xy = slice_gas(**common_arguments_xy, project="mass_weighted_vy")
vx_map_xy = mass_weighted_vx_map_xy / mass_map_xy
vy_map_xy = mass_weighted_vy_map_xy / mass_map_xy

vx_map_xy = vx_map_xy.value
vy_map_xy = vy_map_xy.value
vplanenorm_map_xy = np.sqrt(vx_map_xy ** 2 + vy_map_xy ** 2)
vx_map_xy /= vplanenorm_map_xy
vy_map_xy /= vplanenorm_map_xy

dimx = 512
dimy = 512
new_x = np.linspace(0, dimx, dimx)
new_y = np.linspace(0, dimy, dimy)

step = 20
ax[0, 1].quiver(
    new_x[::step],
    new_y[::step],
    np.transpose(vx_map_xy.reshape((dimx, dimy)))[::step, ::step],
    np.transpose(vy_map_xy.reshape((dimx, dimy)))[::step, ::step],
    color="black",
    scale=1.5 / step,
    scale_units="xy",  # Fixes the arrow length in data coordinates
    pivot="middle"
    # density=1.0,
    # linewidth=0.2,
    # arrowsize=0.4,
)


data.gas.mass_weighted_vx = data.gas.masses * v[:, 0]
data.gas.mass_weighted_vz = data.gas.masses * v[:, 2]
mass_weighted_vx_map = slice_gas(**common_arguments, project="mass_weighted_vx")
mass_weighted_vz_map = slice_gas(**common_arguments, project="mass_weighted_vz")
vx_map = mass_weighted_vx_map / mass_map
vz_map = mass_weighted_vz_map / mass_map

vx_map = vx_map.value
vz_map = vz_map.value
vplanenorm_map = np.sqrt(vx_map ** 2 + vz_map ** 2)
vx_map /= vplanenorm_map
vz_map /= vplanenorm_map


dimx = 512
dimy = 512
new_x = np.linspace(0, dimx, dimx)
new_y = np.linspace(0, dimy, dimy)

step = 20
ax[0, 0].quiver(
    new_x[::step],
    new_y[::step],
    np.transpose(vx_map.reshape((dimx, dimy)))[::step, ::step],
    np.transpose(vz_map.reshape((dimx, dimy)))[::step, ::step],
    color="black",
    scale=1.5 / step,
    scale_units="xy",  # Fixes the arrow length in data coordinates
    pivot="middle"
    # density=1.0,
    # linewidth=0.2,
    # arrowsize=0.4,
)

data.gas.mass_weighted_Bx = data.gas.masses * B[:, 0]
data.gas.mass_weighted_By = data.gas.masses * B[:, 1]
mass_weighted_Bx_map_xy = slice_gas(**common_arguments_xy, project="mass_weighted_Bx")
mass_weighted_By_map_xy = slice_gas(**common_arguments_xy, project="mass_weighted_By")
Bx_map_xy = mass_weighted_Bx_map_xy / mass_map_xy
By_map_xy = mass_weighted_By_map_xy / mass_map_xy
# Bx_map_xy.convert_to_units(1e-7*unyt.g / (unyt.statA * unyt.s * unyt.s)).values
# By_map_xy.convert_to_units(1e-7*unyt.g / (unyt.statA * unyt.s * unyt.s)).values

Bx_map_xy = Bx_map_xy.value
By_map_xy = By_map_xy.value
Bplanenorm_map_xy = np.sqrt(Bx_map_xy ** 2 + By_map_xy ** 2)
Bx_map_xy /= Bplanenorm_map_xy
By_map_xy /= Bplanenorm_map_xy

dimx = 512
dimy = 512
new_x = np.linspace(0, dimx, dimx)
new_y = np.linspace(0, dimy, dimy)

step = 20
ax[1, 1].quiver(
    new_x[::step],
    new_y[::step],
    np.transpose(Bx_map_xy.reshape((dimx, dimy)))[::step, ::step],
    np.transpose(By_map_xy.reshape((dimx, dimy)))[::step, ::step],
    color="black",
    scale=1.5 / step,
    scale_units="xy",  # Fixes the arrow length in data coordinates
    pivot="middle"
    # density=1.0,
    # linewidth=0.2,
    # arrowsize=0.4,
)

data.gas.mass_weighted_Bx = data.gas.masses * B[:, 0]
data.gas.mass_weighted_Bz = data.gas.masses * B[:, 2]
mass_weighted_Bx_map = slice_gas(**common_arguments, project="mass_weighted_Bx")
mass_weighted_Bz_map = slice_gas(**common_arguments, project="mass_weighted_Bz")
Bx_map = mass_weighted_Bx_map / mass_map
Bz_map = mass_weighted_Bz_map / mass_map

Bx_map = Bx_map.value
Bz_map = Bz_map.value
Bplanenorm_map = np.sqrt(Bx_map ** 2 + Bz_map ** 2)
Bx_map /= Bplanenorm_map
Bz_map /= Bplanenorm_map


dimx = 512
dimy = 512
new_x = np.linspace(0, dimx, dimx)
new_y = np.linspace(0, dimy, dimy)

step = 20
ax[1, 0].quiver(
    new_x[::step],
    new_y[::step],
    np.transpose(Bx_map.reshape((dimx, dimy)))[::step, ::step],
    np.transpose(Bz_map.reshape((dimx, dimy)))[::step, ::step],
    color="black",
    scale=1.5 / step,
    scale_units="xy",  # Fixes the arrow length in data coordinates
    pivot="middle"
    # density=1.0,
    # linewidth=0.2,
    # arrowsize=0.4,
)


index = np.argmax(normv)
vpart = normv[index]
rpart = r[index]
Bpart = normB[index]
Berr = np.maximum(h * abs(divB) / normB, 1e-6)
Berrpart = Berr[index]
hpart = h[index]
vpart.convert_to_units(unyt.km / unyt.s)
rpart.convert_to_units(3.08e18 * 1e3 * unyt.cm)
hpart.convert_to_units(3.08e18 * 1e3 * unyt.cm)
Bpart.convert_to_units(1e-7 * unyt.g / (unyt.statA * unyt.s * unyt.s))
print(
    f"Particle {index}, \nwith maximal velocity {vpart.value} km/s, \nmagnetic field {Bpart.value} muG, \nwith error level {Berrpart.value}, \nwith smoothing length {hpart.value} kPc, \nlocated at {rpart.value} kPc"
)

part_pixel_coords = 512 / 2 * (1 + rpart.value / Lslice_kPc)

# ax[1,0].scatter(part_pixel_coords[0],part_pixel_coords[2],color='red',marker='x')


# add panel with infromation about the run
text_common_args = dict(
    fontsize=10, ha="center", va="center", transform=ax[3, 1].transAxes
)

ax[3, 1].text(
    0.5,
    0.8,
    "Cooling halo with spin at time $t=%.2f$ Gyr" % data.metadata.time,
    **text_common_args,
)
ax[3, 1].text(0.5, 0.7, "swift %s" % git.decode("utf-8"), **text_common_args)
ax[3, 1].text(0.5, 0.6, "Branch %s" % gitBranch.decode("utf-8"), **text_common_args)
ax[3, 1].text(0.5, 0.5, hydroScheme.decode("utf-8"), **text_common_args)
ax[3, 1].text(
    0.5,
    0.4,
    kernel.decode("utf-8") + " with $%.2f$ neighbours" % (neighbours),
    **text_common_args,
)
ax[3, 1].text(
    0.5, 0.3, "Artificial diffusion: $%.2f$ " % (artDiffusion), **text_common_args
)
ax[3, 1].text(
    0.5,
    0.2,
    "Dedner Hyp, Hyp_div(v), Par: $%.2f,%.2f,%.2f$ " % (dedHyp, dedHypDivv, dedPar),
    **text_common_args,
)
ax[3, 1].text(
    0.5, 0.1, "Physical resistivity $\eta$: $%.2f$ " % (eta), **text_common_args
)
ax[3, 1].text(0.5, 0.0, "Number of particles $N_p$: $%.0f$ " % (Np), **text_common_args)
ax[3, 1].axis("off")

fig.tight_layout()

plt.savefig(sys.argv[2], dpi=220)
