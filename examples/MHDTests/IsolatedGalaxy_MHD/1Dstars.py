import numpy as np
import h5py
import argparse
import unyt
import matplotlib.pyplot as plt

from swiftsimio import load
from swiftsimio.visualisation.slice import slice_gas
from scipy.stats import binned_statistic

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
n_gas = data.metadata.n_stars

# Retrieve particle attributes of interest
m = data.stars.masses
coords = data.stars.coordinates
coords_center = coords - data.metadata.boxsize / 2
radius = np.linalg.norm(coords_center, axis=1)
velocities = data.stars.velocities

# calculate specific angular momentum j(r)
rotation_axis = np.array([0, 0, 1])
e_phi = np.cross(rotation_axis[None, :], coords_center)
e_phi = e_phi / np.linalg.norm(e_phi, axis=1)[:, None]
v_phi = (
    velocities[:, 0] * e_phi[:, 0]
    + velocities[:, 1] * e_phi[:, 1]
    + velocities[:, 2] * e_phi[:, 2]
)

r_units = 1e3 * PARSEC_IN_CGS * unyt.cm
v_units = 1e5 * unyt.cm / unyt.s

v_phi = v_phi.to(v_units).value

radius = np.linalg.norm(coords_center, axis=1)
radius = radius.to(r_units).value

# Plot maps
plt.rcParams.update({"font.size": 16})

nx = 1
ny = 2
fig, axs = plt.subplots(ny, nx, figsize=((10 * nx, 5 * ny)))
fig.subplots_adjust(hspace=0.1)

# plot velocity
axs[0].scatter(radius, v_phi, s=0.1, color="black", label="particles")
# axs[0].set_yscale("log")

# plot mean, and 1 sigma interval
num_bins = 100
bins = np.linspace(np.min(radius), np.max(radius), num_bins + 1)
bin_centers = 0.5 * (bins[1:] + bins[:-1])

# Compute binned statistics
mean_vals, _, _ = binned_statistic(radius, v_phi, statistic="mean", bins=bins)
std_vals, _, _ = binned_statistic(radius, v_phi, statistic="std", bins=bins)

# Plot binned mean and ±2σ
axs[0].plot(bin_centers, mean_vals, color="red", label="Mean $v_\\phi$")
axs[0].plot(
    bin_centers, mean_vals + 2 * std_vals, color="red", linestyle=":", label="+2σ"
)
axs[0].plot(
    bin_centers, mean_vals - 2 * std_vals, color="red", linestyle=":", label="-2σ"
)

# limits
axs[0].set_ylim(0, 350)
axs[0].set_yticks([50, 100, 150, 200, 250, 300, 350])

axs[0].set_xlim(1e-2, 2e1)

# add labels
axs[0].set_ylabel(r"$v_{\phi}(r)$ $[\rm km \cdot s^{-1}]$")
axs[0].set_xlabel(r"$r$ $[kpc]$")

# add legend
axs[0].legend()

# vertical lines
axs[0].axvline(x=(r_200_cgs * unyt.cm).to(r_units), color="gray", label=r"$R_{200}$")

Ninfo = 1
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
