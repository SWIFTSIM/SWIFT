################################################################################
# This file is part of SWIFT.
# Copyright (c) 2022 Yves Revaz (yves.revaz@epfl.ch)
#               2026 Darwin Roduit (darwin.roduit@epfl.ch)
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published
# by the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

import json
import os

import h5py
import numpy as np
import argparse
from astropy import units
from astropy import constants


def parse_options():

    usage = "usage: %prog [options] file"
    parser = argparse.ArgumentParser(description=usage)

    parser.add_argument(
        "--rho",
        action="store",
        dest="rho",
        type=float,
        default=5,
        help="Cold-phase gas density in atom/cm3 (also the uniform variant's density)",
    )

    parser.add_argument(
        "--mass",
        action="store",
        dest="mass",
        type=float,
        default=0.1,
        help="Cold-phase gas particle mass in solar mass",
    )

    parser.add_argument(
        "--level",
        action="store",
        dest="level",
        type=int,
        default=6,
        help="Resolution level: N = (2**l)**3",
    )

    parser.add_argument(
        "--boxsize",
        action="store",
        dest="boxsize",
        type=float,
        default=None,
        help="Boxzise in kpc",
    )

    parser.add_argument(
        "--star_mass",
        action="store",
        type=float,
        default=29.7,
        help="Mass of each star in M_sun",
    )

    parser.add_argument(
        "--star_type",
        action="store",
        type=str,
        default="single_star",
        choices=["single_star", "continuous_IMF", "SSP"],
        help="Type of the star for the GEAR model.",
    )

    parser.add_argument(
        "-o",
        action="store",
        dest="outputfilename",
        type=str,
        default="box.hdf5",
        help="output filename",
    )

    parser.add_argument(
        "--bulk-temperature-K",
        type=float,
        default=500.0,
        help="Cold-phase gas temperature (K); the uniform variant's only "
        "temperature. 500 K (not the family's usual 1000 K default) keeps "
        "T_hot = T_cold * 4**bin_delta under the 1e4 K mu-switch at "
        "bin_delta=2 (see README).",
    )

    parser.add_argument(
        "--variant",
        choices=["uniform", "twophase"],
        default="twophase",
        help="uniform: single-bin control, every particle at the cold "
        "phase's density/mass/temperature. twophase: the multi-bin test, "
        "split at x=L/2 (see README).",
    )

    parser.add_argument(
        "--bin-delta",
        type=int,
        default=2,
        help="Timestep-bin separation the two phases are built for: "
        "T_hot = T_cold * 4**bin_delta, n_hot = n_cold / 4**bin_delta, "
        "m_hot = m_cold / 4**bin_delta. 2 is the maximum meaningful value: "
        "--limiter clamps neighbour-pair bin differences at 2 "
        "(time_bin_neighbour_max_delta_bin), so a larger split is never "
        "realized on a neighbour pair.",
    )

    parser.add_argument(
        "--allow-clamped-bins",
        action="store_true",
        help="Required to pass --bin-delta > 2: acknowledges the request "
        "will be silently clamped to 2 by the limiter and is not testing "
        "what its own value implies.",
    )

    parser.add_argument(
        "--interface-offset-h",
        type=float,
        default=1.0,
        help="Star A's distance inside the cold phase from the phase "
        "plane at x=L/2, in units of the glass's own median h.",
    )

    parser.add_argument(
        "--v-over-c-hyp",
        type=float,
        default=0.25,
        help="Star velocity as a fraction of the cold phase's own c_hyp, "
        "used only when --star-velocity-km-s is negative (derive).",
    )

    parser.add_argument(
        "--star-velocity-km-s",
        type=float,
        default=-1.0,
        help="Star velocity (+y, all stars) in km/s. Negative (default): "
        "derive from --v-over-c-hyp and the cold phase's own analytic "
        "c_hyp.",
    )

    parser.add_argument(
        "--c-hyp-margin",
        type=float,
        default=0.5,
        help="GEARFeedback:LW_FUV_c_hyp_margin, used only to derive "
        "--star-velocity-km-s (must match the run's own params.yml value).",
    )

    parser.add_argument(
        "--cfl",
        type=float,
        default=0.1,
        help="SPH:CFL_condition, used only to derive --star-velocity-km-s "
        "(must match the run's own params.yml value).",
    )

    parser.add_argument(
        "--n-stars",
        type=int,
        default=3,
        choices=[1, 3],
        help="3 (default, shipped): stars B, A, C of the design (see "
        "README). 1 (debugging only): star A alone.",
    )

    options = parser.parse_args()
    return options


def internal_energy_from_temperature_cgs(T_K, hydrogen_mass_fraction=0.752):
    """u = kB*T / ((gamma-1)*mu*m_p), cgs (erg/g); mirrors src/hydro_properties.c's own formula (gamma=5/3, mu switches at the 1e4 K default ionization threshold)."""
    kB_cgs = 1.380649e-16
    mp_cgs = 1.67262192369e-24
    gamma_minus_one_inv = 1.5  # gamma = 5/3
    H_ionization_temperature_K = 1e4
    if T_K > H_ionization_temperature_K:
        mu = 4.0 / (8.0 - 5.0 * (1.0 - hydrogen_mass_fraction))
    else:
        mu = 4.0 / (1.0 + 3.0 * hydrogen_mass_fraction)
    return kB_cgs * T_K * gamma_minus_one_inv / (mu * mp_cgs)


def sound_speed_km_s_from_temperature(T_K):
    """c_s = sqrt(gamma*(gamma-1)*u), gamma = 5/3; converted from the cgs
    internal energy above to km/s."""
    u_cgs = internal_energy_from_temperature_cgs(T_K)
    gamma = 5.0 / 3.0
    cs_cgs = np.sqrt(gamma * (gamma - 1.0) * u_cgs)  # cm/s
    return cs_cgs / 1e5  # km/s


########################################
# main
########################################

opt = parse_options()

if opt.bin_delta > 2 and not opt.allow_clamped_bins:
    raise RuntimeError(
        f"--bin-delta {opt.bin_delta} > 2 requested without "
        "--allow-clamped-bins: --limiter clamps neighbour-pair timestep-bin "
        "differences at 2 (time_bin_neighbour_max_delta_bin), so a larger "
        "split is never realized on a neighbour pair and this run would not "
        "test what its own --bin-delta implies. Pass --allow-clamped-bins "
        "to run it anyway."
    )

T_cold = opt.bulk_temperature_K
T_hot = T_cold * 4.0**opt.bin_delta
if T_hot >= 1e4:
    raise RuntimeError(
        f"T_hot = {T_hot:.1f} K >= 1e4 K: internal_energy_from_temperature_cgs "
        "switches its mu at the 1e4 K ionization threshold, which would make "
        "the hot/cold sound-speed ratio a non-power-of-two and the realized "
        "bin split non-deterministic (constraint C6). Lower "
        "--bulk-temperature-K or --bin-delta."
    )

# define standard units
UnitMass_in_cgs = 1.988409870698051e43  # 10^10 M_sun in grams
UnitLength_in_cgs = 3.0856775814913673e21  # kpc in centimeters
UnitVelocity_in_cgs = 1e5  # km/s in centimeters per second
UnitCurrent_in_cgs = 1  # Amperes
UnitTemp_in_cgs = 1  # Kelvin
UnitTime_in_cgs = UnitLength_in_cgs / UnitVelocity_in_cgs

UnitMass = UnitMass_in_cgs * units.g
UnitLength = UnitLength_in_cgs * units.cm
UnitTime = UnitTime_in_cgs * units.s
UnitVelocity = UnitVelocity_in_cgs * units.cm / units.s

np.random.seed(1)

# Number of particles
if opt.level is not None:
    N = (2**opt.level) ** 3

# Mean (cold-phase) density, used to size the box exactly as the family
# does: the box is sized from the cold-phase density/mass as if it were
# uniform, then the twophase variant relabels the x>=L/2 half hot.
rho = opt.rho  # atom/cc
rho = rho * constants.m_p / units.cm**3

# Cold-phase gas particle mass
m = opt.mass  # in solar mass
m = m * units.Msun

# Size of the box
if opt.level is not None:
    if opt.boxsize is None:
        M = N * m
        L = (M / rho) ** (1 / 3.0)
    else:
        L = opt.boxsize * units.kpc
        M = rho * L**3
        print(M.to(units.Msun))
        N = int(np.ceil(M.to(units.Msun) / m))
        print(N)
else:
    if opt.boxsize is not None:
        L = opt.boxsize * units.kpc
        M = (rho * L**3).to(units.Msun)
        N = int(np.ceil(M / m))
    else:
        raise RuntimeError(
            "If the opt.level and opt.boxsize are None, then we cannot determine the number of particles"
        )

print("Boxsize                               : {}".format(L.to(units.kpc)))
print("Total mass (cold-phase sizing)        : {}".format(M.to(units.Msun)))
print("Number of particles                   : {}".format(N))

# Convert to code units
m_cold = m.to(UnitMass).value
L = L.to(UnitLength).value
rho_cold = rho.to(UnitMass / UnitLength**3).value

scale = 4.0**opt.bin_delta
m_hot = m_cold / scale
rho_hot = rho_cold / scale

# Pre-relaxed glass (getGlass.sh): random placement's shot noise would
# swamp the hydrostatic two-phase structure this example depends on.
if opt.level is not None and opt.boxsize is None:
    glass_n = 2**opt.level
    glass_filename = f"glassCube_{glass_n}.hdf5"
    if not os.path.exists(glass_filename):
        raise RuntimeError(
            f"Missing {glass_filename}. Run ./getGlass.sh {glass_n} first."
        )
    with h5py.File(glass_filename, "r") as glass:
        pos = glass["/PartType0/Coordinates"][:, :]
        h_glass = glass["/PartType0/SmoothingLength"][:]
    if pos.shape[0] != N:
        raise RuntimeError(
            f"{glass_filename} has {pos.shape[0]} particles, expected {N} "
            f"for level={opt.level}."
        )
    eps = 1e-6
    pos = (pos - pos.min()) / (pos.max() - pos.min() + eps) * L
    h = h_glass * L
else:
    print(
        "No matching glass file for this --boxsize/--level combination; "
        "falling back to a random, non-thermally-relaxed placement."
    )
    pos = np.random.random([N, 3]) * np.array([L, L, L])
    h = np.ones(N) * 3 * L / N ** (1 / 3.0)

h_median = float(np.median(h))
print("Inter-particle distance (code unit)   : {}".format(L / N ** (1 / 3.0)))
print("Median smoothing length (code unit)   : {:.6e}".format(h_median))

vel = np.zeros([N, 3])
ids = np.arange(N)

# Phase assignment (constraint C10: recoverable from Masses alone, since
# TimeBin is not a snapshot field).
if opt.variant == "twophase":
    is_hot = pos[:, 0] >= (L / 2.0)
else:
    is_hot = np.zeros(N, dtype=bool)

mass = np.where(is_hot, m_hot, m_cold).astype(np.float64)
rho_arr = np.where(is_hot, rho_hot, rho_cold).astype(np.float64)

UnitVelocity_cgs2 = UnitVelocity_in_cgs**2
u_cold = internal_energy_from_temperature_cgs(T_cold) / UnitVelocity_cgs2
u_hot = internal_energy_from_temperature_cgs(T_hot) / UnitVelocity_cgs2
u = np.where(is_hot, u_hot, u_cold).astype(np.float64)

n_hot = int(np.sum(is_hot))
n_cold = N - n_hot
print(
    f"Phase populations: cold={n_cold} ({n_cold/N:.1%})  "
    f"hot={n_hot} ({n_hot/N:.1%})"
)

#####################
# Kinematics: derive (or take) the star velocity from the cold phase's own
# analytic c_hyp.
#####################
cs_cold_km_s = sound_speed_km_s_from_temperature(T_cold)
cs_hot_km_s = sound_speed_km_s_from_temperature(T_hot)
kernel_gamma = 1.936492  # Wendland-C2, 3D

# dt = 2*kernel_gamma*CFL*h/(2*c_s); h in kpc, c_s in km/s, UnitTime =
# kpc/(km/s), so this ratio is already in internal time units. dt scales
# with the sound-speed ratio (sqrt(scale) = 4 at bin_delta=2), not with
# `scale` itself (16, the mass/density ratio): computed directly from
# cs_hot rather than dividing dt_cold by `scale`.
dt_cold_analytic = kernel_gamma * opt.cfl * h_median / cs_cold_km_s
dt_hot_analytic = kernel_gamma * opt.cfl * h_median / cs_hot_km_s

c_hyp_cold_analytic = opt.c_hyp_margin * h_median / dt_cold_analytic
c_hyp_hot_analytic = opt.c_hyp_margin * h_median / dt_hot_analytic

if opt.star_velocity_km_s >= 0.0:
    v_star = opt.star_velocity_km_s
else:
    v_star = opt.v_over_c_hyp * c_hyp_cold_analytic

print(f"Sound speed, cold phase (km/s)         : {cs_cold_km_s:.6f}")
print(f"Sound speed, hot phase (km/s)           : {cs_hot_km_s:.6f}")
print(f"dt_cold (analytic, internal)            : {dt_cold_analytic:.6e}")
print(f"dt_hot (analytic, internal)             : {dt_hot_analytic:.6e}")
print(f"c_hyp_cold (analytic, km/s)             : {c_hyp_cold_analytic:.6f}")
print(f"c_hyp_hot (analytic, km/s)              : {c_hyp_hot_analytic:.6f}")
print(f"Star velocity v_star (km/s)             : {v_star:.6f}")
print(f"v_star / c_hyp_cold                     : {v_star/c_hyp_cold_analytic:.4f}")

#####################
# Stars
#####################
G = constants.G
M_star = (opt.star_mass * units.M_sun).to(UnitMass).value

x_B = L / 4.0
x_A = L / 2.0 - opt.interface_offset_h * h_median
x_C = 3.0 * L / 4.0
y_mid = L / 2.0
z_mid = L / 2.0

if opt.n_stars == 3:
    star_names = ["B", "A", "C"]
    star_x = [x_B, x_A, x_C]
else:
    star_names = ["A"]
    star_x = [x_A]

N_star = len(star_x)
pos_star = np.array([[x, y_mid, z_mid] for x in star_x])
vel_star = np.zeros([N_star, 3])
vel_star[:, 1] = v_star

for name, x in zip(star_names, star_x):
    dist_to_plane_h = (x - L / 2.0) / h_median
    print(
        f"Star {name}: x={x:.6e} ({x/h_median:.3f} h)  "
        f"distance to phase plane = {dist_to_plane_h:.3f} h"
    )

M_star_arr = np.full(N_star, M_star)
h_star = np.ones(N_star) * 3 * L / N ** (1 / 3.0)  # Same as the gas particles

if opt.star_type == "single_star":
    star_type = 0
elif opt.star_type == "continuous_IMF":
    star_type = 1
elif opt.star_type == "SSP":
    star_type = 2
else:
    raise ValueError(f"Star type {opt.star_type} is not known!")

star_particle_type = np.ones(N_star) * star_type
star_id = N + np.arange(1, N_star + 1)
star_birth_time = np.zeros(N_star)

#####################
# Finally write the ICs in the file
#####################

fileOutput = h5py.File(opt.outputfilename, "w")
print("{} saved.".format(opt.outputfilename))

# Header
grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [L, L, L]
grp.attrs["NumPart_Total"] = [N, 0, 0, 0, N_star, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, N_star, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFileOutputsPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3

# Units
grp = fileOutput.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = UnitLength_in_cgs
grp.attrs["Unit mass in cgs (U_M)"] = UnitMass_in_cgs
grp.attrs["Unit time in cgs (U_t)"] = UnitTime_in_cgs
grp.attrs["Unit current in cgs (U_I)"] = UnitCurrent_in_cgs
grp.attrs["Unit temperature in cgs (U_T)"] = UnitTemp_in_cgs

# Write Gas particle group. FUVSpecificEnergy/LWSpecificEnergy are not
# seeded: the field starts at zero and is built by the stars.
grp = fileOutput.create_group("/PartType0")
grp.create_dataset("Coordinates", data=pos, dtype="d")
grp.create_dataset("Velocities", data=vel, dtype="f")
grp.create_dataset("Masses", data=mass, dtype="f")
grp.create_dataset("SmoothingLength", data=h, dtype="f")
grp.create_dataset("InternalEnergy", data=u, dtype="f")
grp.create_dataset("ParticleIDs", data=ids, dtype="L")
grp.create_dataset("Densities", data=rho_arr, dtype="f")

# Write star particle group
grp = fileOutput.create_group("/PartType4")
grp.create_dataset("Coordinates", data=pos_star, dtype="d")
grp.create_dataset("Velocities", data=vel_star, dtype="f")
grp.create_dataset("Masses", data=M_star_arr, dtype="f")
grp.create_dataset("ParticleIDs", data=star_id, dtype="L")
grp.create_dataset("SmoothingLength", data=h_star, dtype="f")
grp.create_dataset("BirthMass", data=M_star_arr, dtype="f")
grp.create_dataset("BirthTime", data=star_birth_time, dtype="f")
grp.create_dataset("StellarParticleType", data=star_particle_type, dtype="i")

fileOutput.close()

# Sidecar: the check script reads the derived quantities from here instead
# of re-deriving them.
sidecar = dict(
    L=L,
    h_median=h_median,
    m_cold=m_cold,
    m_hot=m_hot,
    T_cold=T_cold,
    T_hot=T_hot,
    bin_delta=opt.bin_delta,
    v_star=v_star,
    dt_cold_analytic=dt_cold_analytic,
    dt_hot_analytic=dt_hot_analytic,
    c_hyp_cold_analytic=c_hyp_cold_analytic,
    c_hyp_hot_analytic=c_hyp_hot_analytic,
    star_positions={n: [x, y_mid, z_mid] for n, x in zip(star_names, star_x)},
    star_names=star_names,
    variant=opt.variant,
    cs_cold_km_s=cs_cold_km_s,
    cs_hot_km_s=cs_hot_km_s,
    interface_offset_h=opt.interface_offset_h,
)
with open("multibin_ic.json", "w") as f:
    json.dump(sidecar, f, indent=2)
print("multibin_ic.json saved.")
