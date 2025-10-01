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

radius = np.linalg.norm(coords_center, axis=1)

# Plot maps
plt.rcParams.update({"font.size": 16})

nx = 1
ny = 4
fig, axs = plt.subplots(ny, nx, figsize=((10 * nx, 5 * ny)))
fig.subplots_adjust(hspace=0.1)

# plot density
axs[0].scatter(radius.to(r_units), n_H.to(nH_units), s=0.1, color="black")
axs[0].set_yscale("log")
axs[0].set_yticks(np.logspace(-7, 2, 10))
axs[0].set_ylim(1e-7, 1e2)

# plot Pressure
axs[1].scatter(radius.to(r_units), P.to(P_units), s=0.1, color="black")
axs[1].set_yscale("log")

# plot specific angular momentum
axs[2].scatter(radius.to(r_units), j.to(j_units), s=0.1, color="black")
axs[2].set_yscale("log")

# NFW-like gas density profile
def rho_r(r_value, f_b, M_200_cgs, r_200_cgs, c_200):
    rho_0 = (
        M_200_cgs
        / (np.log(1 + c_200) - c_200 / (1 + c_200))
        / (4 * np.pi * r_200_cgs ** 3 / c_200 ** 3)
    )
    result_cgs = rho_0 * f_b / (c_200 * r_value * (1 + c_200 * r_value) ** 2)
    # Apply density cut
    rho_max_cgs = nH_max_cgs * m_H_cgs
    result_cgs = np.array(result_cgs)
    result_cgs[result_cgs > rho_max_cgs] = rho_max_cgs
    return result_cgs


# NFW-like gas mass inside a sphere with radius R
def Mgas_r(r_value, f_b, M_200_cgs, r_200_cgs, c_200):
    M_0 = M_200_cgs / (np.log(1 + c_200) - c_200 / (1 + c_200))
    return (
        M_0
        * f_b
        * (np.log(1 + c_200 * r_value) - c_200 * r_value / (1 + c_200 * r_value))
    )


# NFW Gravitational acceleration
def a_NFW(r_value, M_200_cgs, r_200_cgs, c_200):
    a_pref = (
        CONST_G_CGS
        * M_200_cgs
        / (np.log(1 + c_200) - c_200 / (1 + c_200))
        / r_200_cgs ** 2
    )
    return (
        a_pref
        * ((r_value / (r_value + 1 / c_200)) - np.log(1 + c_200 * r_value))
        / r_value ** 2
    )


# Integrate rho_gas*a_NFW
def integrate(r_min, r_max, f_b, M_200_cgs, r_200_cgs, c_200, Nsteps=10000):
    # Perform the integration
    r_range = np.linspace(r_min, r_max, Nsteps)
    dr = np.abs((r_max - r_min) / Nsteps)
    integrands = rho_r(r_range, f_b, M_200_cgs, r_200_cgs, c_200) * a_NFW(
        r_range, M_200_cgs, r_200_cgs, c_200
    )
    result_cgs = np.sum(integrands * dr) * r_200_cgs
    return result_cgs


# NFW-like gas hydrostatic equilibrium internal energy profile
def u_vs_r(P0_cgs, r_value, r_max, f_b, M_200_cgs, r_200_cgs, c_200):
    result_cgs = (
        (P0_cgs - integrate(r_value, r_max, f_b, M_200_cgs, r_200_cgs, c_200))
        / (gamma - 1)
        / rho_r(r_value, f_b, M_200_cgs, r_200_cgs, c_200)
    )
    return result_cgs


# NFW-like gas hydrostatic equilibrium pressure profile
def P_vs_r(P0_cgs, r_value, r_max, f_b, M_200_cgs, r_200_cgs, c_200):
    result_cgs = P0_cgs - integrate(r_value, r_max, f_b, M_200_cgs, r_200_cgs, c_200)
    return result_cgs


# NFW-like gas hydrostatic equilibrium temperature profile
def T_vs_r(P0_cgs, r_value, r_max, f_b, M_200_cgs, r_200_cgs, c_200):
    result_cgs = (
        u_vs_r(P0_cgs, r_value, r_max, f_b, M_200_cgs, r_200_cgs, c_200)
        * (gamma - 1)
        * m_H_cgs
        / kb_cgs
    )
    return result_cgs


# define specific angular momentum distribution
def j(r_value, j_max, s, f_b, M_200_cgs, r_200_cgs, c_200):
    return (
        j_max
        * (Mgas_r(r_value, f_b, M_200_cgs, r_200_cgs, c_200) / (M_200_cgs * f_b)) ** s
    )


def jsimple(r_value, j_max, s, f_b, M_200_cgs, r_200_cgs, c_200):
    return j_max * r_value ** s


rmax = (r_200_cgs * unyt.cm).to(radius.units)  # np.max(radius)
rmin = np.min(radius)
r_analytic = (
    np.logspace(np.log10(rmin.value), np.log10(rmax.value), 1000) * radius.units
)
# r_analytic = np.linspace(rmin.value,rmax.value,10000)*radius.units
rho_analytic = rho_r(
    (r_analytic / (r_200_cgs * unyt.cm)).value, f_b, M_200_cgs, r_200_cgs, c_200
)
Mgas_analytic = Mgas_r(
    (r_analytic / (r_200_cgs * unyt.cm)).value, f_b, M_200_cgs, r_200_cgs, c_200
)
nH_analytic = rho_analytic / m_H_cgs * nH_units


T0_cgs = (
    T_200_cgs
)  # 1e5 # gas temperature on the edge of the box (if we want to set this manually)
rho0_cgs = rho_r((rmax / (r_200_cgs * unyt.cm)).value, f_b, M_200_cgs, r_200_cgs, c_200)

P0_cgs = rho0_cgs * kb_cgs * T0_cgs / m_H_cgs  # gas pressure on the edge of the box
P_analytic = np.array(
    [
        P_vs_r(
            P0_cgs,
            (r_analytic[i] / (r_200_cgs * unyt.cm)).value,
            (rmax / (r_200_cgs * unyt.cm)).value,
            f_b,
            M_200_cgs,
            r_200_cgs,
            c_200,
        )
        for i in range(len(r_analytic))
    ]
) * (
    unyt.g / (unyt.cm * unyt.s ** 2)
)  # gas particle internal energies

# Normalize to Bullock
mask_r200 = radius <= (r_200_cgs * unyt.cm)
Jtot = np.sum(omega[mask_r200] * axis_dist[mask_r200] ** 2 * m[mask_r200])
jsp = Jtot / np.sum(m[mask_r200])

jsp_B = (
    (np.sqrt(2) * v_200_cgs * (unyt.cm / unyt.s)) * r_200_cgs * (unyt.cm) * spin_lambda
)
print("Total spin ratio: ", jsp_B.to(j_units).value / jsp.to(j_units).value)

############

j_1_analytic = j(
    (r_analytic / (r_200_cgs * unyt.cm)).value, 1, 1, f_b, M_200_cgs, r_200_cgs, c_200
)
j_1_analytic_simple = jsimple(
    (r_analytic / (r_200_cgs * unyt.cm)).value, 1, 2, f_b, M_200_cgs, r_200_cgs, c_200
)

int_dOmega = (
    8 * np.pi / 3
)  # integral over angles from axial distance squared: int Sin^2 (theta)dphi d Cos (theta)

jsp_1_analytic = (
    np.sum(
        j_1_analytic[:-1]
        * rho_analytic[:-1]
        * (unyt.g / unyt.cm ** 3)
        * int_dOmega
        * r_analytic[:-1] ** 2
        * np.diff(r_analytic)
    )
) / (
    np.sum(
        rho_analytic[:-1]
        * (unyt.g / unyt.cm ** 3)
        * 4
        * np.pi
        * r_analytic[:-1] ** 2
        * np.diff(r_analytic)
    )
)
jsp_1_analytic_simple = (
    np.sum(
        j_1_analytic_simple[:-1]
        * rho_analytic[:-1]
        * (unyt.g / unyt.cm ** 3)
        * int_dOmega
        * r_analytic[:-1] ** 2
        * np.diff(r_analytic)
    )
) / (
    np.sum(
        rho_analytic[:-1]
        * (unyt.g / unyt.cm ** 3)
        * 4
        * np.pi
        * r_analytic[:-1] ** 2
        * np.diff(r_analytic)
    )
)

j_analytic = j_1_analytic / jsp_1_analytic * jsp_B

j_analytic_simple = j_1_analytic_simple / jsp_1_analytic_simple * jsp_B


#
# plot density
axs[0].plot(
    r_analytic.to(r_units),
    nH_analytic.to(nH_units),
    color="red",
    label=r"$\rho_{\rm NFW}^{gas}(r)$",
)
axs[0].set_yscale("log")
axs[0].set_xscale("log")
axs[0].set_xlim([rmin, rmax])
axs[0].legend()

# plot Pressure
axs[1].plot(
    r_analytic.to(r_units),
    P_analytic.to(P_units),
    color="red",
    label=r"$P_{hyd. eq.}(r)$",
)
axs[1].set_yscale("log")
axs[1].set_xscale("log")
axs[1].set_xlim([rmin, rmax])
axs[1].legend()

# plot j(r)
axs[2].plot(
    r_analytic.to(r_units),
    j_analytic.to(j_units),
    color="red",
    label=r"$\rm j(r)\sim M(r)$",
)
axs[2].plot(
    r_analytic.to(r_units),
    j_analytic_simple.to(j_units),
    color="red",
    linestyle="dashed",
    label=r"$\rm j(r)\sim r^2$",
)
axs[2].set_xscale("log")
axs[2].set_xlim([rmin, rmax])
axs[2].legend()

# limits
axs[0].set_ylim(1e-7, 1e2)
axs[0].set_yticks(np.logspace(-7, 2, 10))
axs[1].set_ylim(1e1, 1e9)
axs[1].set_yticks(np.logspace(1, 9, 9))
axs[2].set_ylim(1e18, 1e27)
axs[2].set_yticks(np.logspace(18, 27, 10))

axs[0].set_xlim(1e-2, 1e3)
axs[1].set_xlim(1e-2, 1e3)
axs[2].set_xlim(1e-2, 1e3)

# add labels
axs[0].set_ylabel(r"$n_H(r)$ $[cm^{-3}]$")
axs[1].set_ylabel(r"$P(r)$ $[\rm M_{sol} \cdot pc^{-1} \cdot Gyr^{-2}]$")
axs[2].set_ylabel(r"$j(r)$ $[\rm kpc^{2} \cdot Gyr^{-1} ]$")

# vertical lines
axs[0].axvline(x=(r_200_cgs * unyt.cm).to(r_units), color="gray", label=r"$R_{200}$")
axs[1].axvline(x=(r_200_cgs * unyt.cm).to(r_units), color="gray", label=r"$R_{200}$")
axs[2].axvline(x=(r_200_cgs * unyt.cm).to(r_units), color="gray", label=r"$R_{200}$")

Ninfo = 3
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
