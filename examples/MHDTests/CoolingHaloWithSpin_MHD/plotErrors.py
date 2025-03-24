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
mu0 = 1.25663706127e-1 * unyt.g * unyt.cm / (unyt.s ** 2 * unyt.A ** 2)

# First create a mass-weighted temperature dataset
r = data.gas.coordinates - center
v = data.gas.velocities
norm_r = np.sqrt(r[:, 0] ** 2 + r[:, 1] ** 2 + r[:, 2] ** 2)
vr = np.sum(v * r, axis=1) / norm_r
normv = np.sqrt(v[:, 0] ** 2 + v[:, 1] ** 2 + v[:, 2] ** 2)
B = data.gas.magnetic_flux_densities 
#B = v.value * 0.0 * 1e-7*unyt.g / (unyt.statA * unyt.s * unyt.s)
#J = data.gas.magnetic_flux_curl
normB = np.sqrt(B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2)
#normJ = np.sqrt(J[:, 0] ** 2 + J[:, 1] ** 2 + J[:, 2] ** 2)
divB = data.gas.magnetic_divergences
#divB = norm_r.value * 0.0 * 1e-7*unyt.g / (unyt.statA * unyt.s * unyt.s* unyt.cm)
h = data.gas.smoothing_lengths
rho = data.gas.densities


R0 = data.gas.r0
R1 = data.gas.r1
R2 = data.gas.r2
OW = data.gas.owtriggers

# Normalize errors
R0[R0<1e-2] = 1e-2
R1[R1<1e-2] = 1e-2
R2[R2<1e-2] = 1e-2
OW[OW<1e-2] = 1e-2



Np = len(h)
print('Number of particles %E' % len(data.gas.masses))
print('Total mass %E Msol' % np.sum(data.gas.masses.to(unyt.Msun).value))
print('Gas particle mass %E Msol' % np.mean(data.gas.masses.to(unyt.Msun).value))

data.gas.mass_weighted_R0 = data.gas.masses * R0
data.gas.mass_weighted_R1 = data.gas.masses * R1
data.gas.mass_weighted_R2 = data.gas.masses * R2
data.gas.mass_weighted_OW = data.gas.masses * OW


"""
rotation_matrix = [[1,0,0],
                   [0,1,0],
                   [0,0,1]]
"""
rotation_matrix = [[0, 0, 1], [0, 1, 0], [-1, 0, 0]]

# set plot area
Lslice_kPc = 20.0 #2*21.5 
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
mass_map = slice_gas(**common_arguments, project="masses")

mass_weighted_R0_map_xy = slice_gas(**common_arguments_xy, project="mass_weighted_R0")
mass_weighted_R1_map_xy = slice_gas(**common_arguments_xy, project="mass_weighted_R1")
mass_weighted_R2_map_xy = slice_gas(**common_arguments_xy, project="mass_weighted_R2")
mass_weighted_OW_map_xy = slice_gas(**common_arguments_xy, project="mass_weighted_OW")
mass_weighted_R0_map = slice_gas(**common_arguments, project="mass_weighted_R0")
mass_weighted_R1_map = slice_gas(**common_arguments, project="mass_weighted_R1")
mass_weighted_R2_map = slice_gas(**common_arguments, project="mass_weighted_R2")
mass_weighted_OW_map = slice_gas(**common_arguments, project="mass_weighted_OW")

#mass_weighted_J_map = slice_gas(**common_arguments, project="mass_weighted_J")
R0_map_xy = mass_weighted_R0_map_xy / mass_map_xy
R1_map_xy = mass_weighted_R1_map_xy / mass_map_xy
R2_map_xy = mass_weighted_R2_map_xy / mass_map_xy
OW_map_xy = mass_weighted_OW_map_xy / mass_map_xy

R0_map = mass_weighted_R0_map / mass_map
R1_map = mass_weighted_R1_map / mass_map
R2_map = mass_weighted_R2_map / mass_map
OW_map = mass_weighted_OW_map / mass_map

nx = 2
ny = 5
fig, ax = plt.subplots(ny, nx, sharey=True, figsize=(5*nx, 5*ny))

#fig.suptitle("Plotting at %.3f free fall times" % toplot)

a00 = ax[0, 0].contourf(
    np.log10(R0_map.value),
    cmap="plasma",#"plasma",#"gist_stern",#"RdBu",
    extend="both",
    levels=np.linspace(-2, 0, 100),
)

a01 = ax[0, 1].contourf(
    np.log10(R0_map_xy.value).T,
    cmap="plasma",#"plasma",#"gist_stern", #"RdBu",
    extend="both",
    levels=np.linspace(-2, 0, 100),
)
a10 = ax[1, 0].contourf(
    np.log10(R1_map.value),
    cmap="plasma",#"viridis", #"nipy_spectral_r",
    extend="both",
    levels=np.linspace(-2, 0, 100),
)
a11 = ax[1, 1].contourf(
    np.log10(R1_map_xy.value).T,
    cmap="plasma",#"viridis", #"nipy_spectral_r",
    extend="both",
    levels=np.linspace(-2, 0, 100),
)

a20 = ax[2, 0].contourf(
    np.log10(R2_map.value),
    cmap="plasma",#"viridis", #"nipy_spectral_r",
    extend="both",
    levels=np.linspace(-2, 0, 100),
)

a21 = ax[2, 1].contourf(
    np.log10(R2_map_xy.value).T,
    cmap="plasma",#"viridis", #"nipy_spectral_r",
    extend="both",
    levels=np.linspace(-2, 0, 100),
)

a30 = ax[3, 0].contourf(
    np.log10(OW_map.value),
    cmap="plasma",#"viridis", #"nipy_spectral_r",
    extend="both",
    levels=np.linspace(-2, 0, 100),
)

a31 = ax[3, 1].contourf(
    np.log10(OW_map_xy.value).T,
    cmap="plasma",#"viridis", #"nipy_spectral_r",
    extend="both",
    levels=np.linspace(-2, 0, 100),
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
        #ax[ii, jj].set_aspect("equal")
        ax[ii,jj].set_aspect('equal', adjustable='box')
        ax[ii,jj].set_xlim(0, 511)
        ax[ii,jj].set_ylim(0, 511)


ticks_error =  [-2,-1,0]
cbar1 = plt.colorbar(
    a00, ax=ax[0, 0], fraction=0.046, pad=0.04, ticks=ticks_error
)
cbar1.set_label(r"$\mathrm{log}_{10} R_0 \; $")
cbar2 = plt.colorbar(
    a01, ax=ax[0, 1], fraction=0.046, pad=0.04, ticks=ticks_error
)
ax[0, 1].set_ylabel(r"$y$ [kPc]")
cbar2.set_label(r"$\mathrm{log}_{10} R_0 \; $")
ax[1, 1].set_ylabel(r"$y$ [kPc]")
cbar3 = plt.colorbar(
    a11, ax=ax[1, 0], fraction=0.046, pad=0.04, ticks=ticks_error
)
cbar3.set_label(r"$\mathrm{log}_{10} R_1 \;$")
cbar4 = plt.colorbar(
    a11, ax=ax[1, 1], fraction=0.046, pad=0.04, ticks=ticks_error
)
cbar4.set_label(r"$\mathrm{log}_{10} R_1 \;$")

ax[2, 1].set_ylabel(r"$y$ [kPc]")

cbar5 = plt.colorbar(a20, ax=ax[2, 0], fraction=0.046, pad=0.04,  ticks=ticks_error)
cbar5.set_label(r"$\mathrm{log}_{10} R_2 \;$")

cbar6 = plt.colorbar(a21, ax=ax[2, 1], fraction=0.046, pad=0.04, ticks=ticks_error)
cbar6.set_label(r"$\mathrm{log}_{10} R_2 \;$")

ticksOW=[-1,0,1]
cbar7 = plt.colorbar(a30, ax=ax[3, 0], fraction=0.046, pad=0.04,ticks=ticksOW)
cbar7.set_label(r"$\mathrm{log}_{10} OW \;$")

cbar8 = plt.colorbar(a31, ax=ax[3, 1], fraction=0.046, pad=0.04,ticks=ticksOW)
cbar8.set_label(r"$\mathrm{log}_{10} OW \;$")

Ninfo = 4

# add panel with infromation about the run
text_common_args = dict(
    fontsize=10, ha="center", va="center", transform=ax[Ninfo, 1].transAxes
)

ax[Ninfo, 1].text(
    0.5,
    0.8,
    "Cooling halo with spin at time $t=%.2f$ Gyr" % data.metadata.time,
    **text_common_args,
)
ax[Ninfo, 1].text(0.5, 0.7, "swift %s" % git.decode("utf-8"), **text_common_args)
ax[Ninfo, 1].text(0.5, 0.6, "Branch %s" % gitBranch.decode("utf-8"), **text_common_args)
ax[Ninfo, 1].text(0.5, 0.5, hydroScheme.decode("utf-8"), **text_common_args)
ax[Ninfo, 1].text(
    0.5,
    0.4,
    kernel.decode("utf-8") + " with $%.2f$ neighbours" % (neighbours),
    **text_common_args,
)
ax[Ninfo, 1].text(
    0.5, 0.3, "Artificial diffusion: $%.2f$ " % (artDiffusion), **text_common_args
)
ax[Ninfo, 1].text(
    0.5,
    0.2,
    "Dedner Hyp, Hyp_div(v), Par: $%.2f,%.2f,%.2f$ " % (dedHyp, dedHypDivv, dedPar),
    **text_common_args,
)
ax[Ninfo, 1].text(
    0.5, 0.1, "Physical resistivity $\eta$: $%.2f$ " % (eta), **text_common_args
)
ax[Ninfo, 1].text(
    0.5, 0.0, "Number of particles $N_p$: $%.0f$ " % (Np), **text_common_args
)
ax[Ninfo, 1].axis("off")

fig.tight_layout()

plt.savefig(sys.argv[2],dpi=220)

