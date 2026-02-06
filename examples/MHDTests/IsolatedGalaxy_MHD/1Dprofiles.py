import numpy as np
import h5py
import argparse
import unyt
import matplotlib.pyplot as plt

from swiftsimio import load
from swiftsimio.visualisation.slice import slice_gas

# Some constants cgs
PARSEC_IN_CGS = 3.0856776e18
KM_PER_SEC_IN_CGS = 1.0e5
CONST_G_CGS = 6.672e-8
MSOL_IN_CGS = 1.9891e33  # Solar mass
kb_cgs = 1.38e-16  # boltzmann constant
m_H_cgs = 1.68e-24  # atomic hydrogen mass
GYR_IN_CGS = 3.1536e16  # gigayear
# First set unit velocity and then the circular velocity parameter for the isothermal potential
const_unit_velocity_in_cgs = 1.0e5  # kms^-1

h_dim = 0.704  # 0.681  # 0.67777  # hubble parameter
H_0_cgs = 100.0 * h_dim * KM_PER_SEC_IN_CGS / (1.0e6 * PARSEC_IN_CGS)

# From this we can find the virial radius, the radius within which the average density of the halo is
# 200. * the mean matter density

# Set M200 and get R200 and V200
f_b = 0.17
c_200 = 7.2
spin_lambda = 0.05
nH_max_cgs = 1e2
M_200_cgs = 1e12 * MSOL_IN_CGS
rhoc_cgs = 3 * H_0_cgs ** 2 / (8 * np.pi * CONST_G_CGS)
r_200_cgs = (3 * M_200_cgs / (4 * np.pi * rhoc_cgs * 200)) ** (1 / 3)
v_200_cgs = np.sqrt(CONST_G_CGS * M_200_cgs / r_200_cgs)
v_200 = v_200_cgs / const_unit_velocity_in_cgs
T_200_cgs = m_H_cgs * v_200_cgs ** 2 / (2 * kb_cgs)

const_unit_mass_in_cgs = M_200_cgs
const_unit_length_in_cgs = r_200_cgs
# Derived quantities
const_unit_time_in_cgs = const_unit_length_in_cgs / const_unit_velocity_in_cgs
const_G = (
    CONST_G_CGS
    * const_unit_mass_in_cgs
    * const_unit_time_in_cgs
    * const_unit_time_in_cgs
    / (const_unit_length_in_cgs * const_unit_length_in_cgs * const_unit_length_in_cgs)
)


# Parse command line arguments
argparser = argparse.ArgumentParser()
argparser.add_argument("input")
argparser.add_argument("output")
args = argparser.parse_args()

# Where we want to slice the slice
y0 = 0.5

# Load snapshot
filename = args.input
data = load(filename)

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
n_gas = data.metadata.n_gas

# Retrieve particle attributes of interest
rho = data.gas.densities
m = data.gas.masses
coords = data.gas.coordinates
coords_center = coords - data.metadata.boxsize / 2
radius = np.linalg.norm(coords_center, axis=1)
velocities = data.gas.velocities

# calculate specific angular momentum j(r)
rotation_axis = np.array([0, 0, 1])
e_phi = np.cross(rotation_axis[None, :], coords_center)
e_phi = e_phi / np.linalg.norm(e_phi, axis=1)[:, None]
v_phi = (
    velocities[:, 0] * e_phi[:, 0]
    + velocities[:, 1] * e_phi[:, 1]
    + velocities[:, 2] * e_phi[:, 2]
)
e_r = np.cross(e_phi, rotation_axis[None, :])
axis_dist = (
    coords_center[:, 0] * e_r[:, 0]
    + coords_center[:, 1] * e_r[:, 1]
    + coords_center[:, 2] * e_r[:, 2]
)
j = v_phi * radius ** 2 / axis_dist
omega = v_phi / axis_dist
Jtot = np.sum(omega * axis_dist ** 2 * m)

P = data.gas.pressures
T = data.gas.temperatures
B = data.gas.magnetic_flux_densities
Bx, By = B[:, 0], B[:, 1]
normB = np.sqrt(B[:, 0] ** 2 + B[:, 1] ** 2 + B[:, 2] ** 2)
divB = data.gas.magnetic_divergences
h = data.gas.smoothing_lengths
errB = np.log10(h * abs(divB) / normB)


rho_units = unyt.g * unyt.cm ** (-3)
r_units = 1e3 * PARSEC_IN_CGS * unyt.cm
P_units = (
    MSOL_IN_CGS / PARSEC_IN_CGS / GYR_IN_CGS ** 2 * (unyt.g / (unyt.cm * unyt.s ** 2))
)
j_units = (1e3 * PARSEC_IN_CGS) / GYR_IN_CGS * (unyt.cm ** 2 / unyt.s)
nH_units = unyt.cm ** (-3)

rho.convert_to_units(rho_units)
P.convert_to_units(P_units)
j.convert_to_units(j_units)
n_H = rho / (m_H_cgs * unyt.g)

# filter galaxy
z = coords_center[:, 2]
r_xy = np.linalg.norm(coords_center[:, :2], axis=1)
mask_galaxy = (r_xy.to(r_units) < 30) & (
    np.abs(z).to(r_units) < 2
)  # r_xy < 10 & |z|< 1 galaxy

print(sum(mask_galaxy))

radius = np.linalg.norm(coords_center, axis=1)

# Plot maps
plt.rcParams.update({"font.size": 16})

nx = 1
ny = 3
fig, axs = plt.subplots(ny, nx, figsize=((10 * nx, 5 * ny)))
fig.subplots_adjust(hspace=0.1)

# plot density
axs[0].scatter(
    radius[~mask_galaxy].to(r_units),
    rho[~mask_galaxy].to(rho_units),
    s=0.1,
    color="black",
)
axs[0].set_yscale("log")
# axs[0].set_yticks(np.logspace(-7, 2, 10))
# axs[0].set_ylim(1e-7, 1e2)

# plot density
axs[1].scatter(radius[~mask_galaxy].to(r_units), T[~mask_galaxy], s=0.1, color="black")
axs[1].set_yscale("log")


def rho(r, pars):
    """ 
    Calculate the density profile at given positions
    param: r - distance to the center
    param: pars - parameter dict for density profile.
        pars['profile'] - name of the profile
        other fields - parameters specific for the profile
    return: the density at the given positions
    """

    # Define density profiles
    if pars["profile"] == "uniform":
        # Constant density profile
        rho0 = pars["rho0"]
        density = rho0 * np.ones_like(r)
    elif pars["profile"] == "isothermal":
        # Isothermal profile
        rho0 = pars["rho0"]
        r0 = pars["r0"]
        density = rho0 / (1 + (r / r0) ** 2)
    elif pars["profile"] == "beta":
        # Beta profile
        rho0 = pars["rho0"]
        r0 = pars["r0"]
        beta = pars["beta"]
        density = rho0 * (1 + (r / r0) ** 2) ** (-3 * beta / 2)

    return density


r_analytic = np.logspace(-3, 2, 1000) * r_units

parameters = []
rho0 = 5e-26 * rho_units
R_max = 50 * r_units  # maximal radius in kpc


pars = {
    "density": {"profile": "beta", "rho0": rho0, "r0": 0.33 * r_units, "beta": 2 / 3}
}
# form parameter dict

rho_analytic = rho(r_analytic, pars["density"])

# print(rho_analytic[12:])

# plot density
axs[0].plot(
    r_analytic.to(r_units),
    rho_analytic.to(rho_units),
    color="red",
    label=r"$\rho_{\rm \beta}(r)$",
)
axs[0].set_yscale("log")
axs[0].set_xscale("log")
# axs[0].set_xlim([rmin, rmax])
axs[0].legend()


axs[1].set_yscale("log")
axs[1].set_xscale("log")
# axs[0].set_xlim([rmin, rmax])
axs[1].legend()


# limits
axs[0].set_ylim(1e-30, 1e-24)
axs[0].set_yticks(np.logspace(-30, -24, 4))
axs[0].set_xlim(1e0, 1e2)

axs[1].set_ylim(1e2, 1e8)
axs[1].set_yticks(np.logspace(2, 8, 7))
axs[1].set_xlim(1e0, 1e2)


# add labels
axs[0].set_ylabel(r"$\rho(r)$ $[g/cm^{-3}]$")

# vertical lines
axs[0].axvline(x=(r_200_cgs * unyt.cm).to(r_units), color="gray", label=r"$R_{200}$")

Ninfo = 2
text_common_args = dict(
    fontsize=10, ha="center", va="center", transform=axs[Ninfo].transAxes
)

axs[Ninfo].text(
    0.5, 0.8, r"Cooling Halo with Spin, IC quality check.", **text_common_args
)
axs[Ninfo].text(0.5, 0.7, r"SWIFT %s" % git.decode("utf-8"), **text_common_args)
axs[Ninfo].text(0.5, 0.6, r"Branch %s" % gitBranch.decode("utf-8"), **text_common_args)
axs[Ninfo].text(0.5, 0.5, hydroScheme.decode("utf-8"), **text_common_args)
axs[Ninfo].text(
    0.5,
    0.4,
    kernel.decode("utf-8") + r" with $%.2f$ neighbours" % (neighbours),
    **text_common_args,
)
axs[Ninfo].text(
    0.5, 0.3, r"Dimenstionless hubble constant: $%.3f$ " % h_dim, **text_common_args
)
axs[Ninfo].text(
    0.5, 0.0, r"Number of particles: $N_p$: $%.0f$ " % (n_gas), **text_common_args
)

axs[Ninfo].axis("off")


for ax in axs:
    ax.minorticks_on()
    ax.grid()

fig.tight_layout()
plt.savefig(args.output)
