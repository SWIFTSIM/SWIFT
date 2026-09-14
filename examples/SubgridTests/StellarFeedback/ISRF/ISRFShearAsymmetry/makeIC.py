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

import os
import h5py
import numpy as np
import argparse
from astropy import units
from astropy import constants


class store_as_array(argparse._StoreAction):
    """Provides numpy array as argparse arguments."""

    def __call__(self, parser, namespace, values, option_string=None):
        values = np.array(values)
        return super().__call__(parser, namespace, values, option_string)


def parse_options():

    usage = "usage: %prog [options] file"
    parser = argparse.ArgumentParser(description=usage)

    parser.add_argument(
        "--rho", action="store", dest="rho", type=float, default=5,
        help="Mean gas density in atom/cm3",
    )
    parser.add_argument(
        "--mass", action="store", dest="mass", type=float, default=0.1,
        help="Gas particle mass in solar mass",
    )
    parser.add_argument(
        "--level", action="store", dest="level", type=int, default=6,
        help="Resolution level: N = (2**l)**3",
    )
    parser.add_argument(
        "--boxsize", action="store", dest="boxsize", type=float, default=None,
        help="Boxzise in kpc",
    )
    parser.add_argument(
        "-o", action="store", dest="outputfilename", type=str, default="box.hdf5",
        help="output filename",
    )

    parser.add_argument(
        "--variant", type=str, default="shear", choices=["shear", "contrast"],
        help="shear: uniform gas, sheared velocity only. contrast: also a "
        "density/temperature contrast across the shear layer.",
    )
    parser.add_argument(
        "--v-shear-km-s", type=float, default=1.0,
        help="Full v_rel between the two streams (km/s).",
    )
    parser.add_argument(
        "--layer-width-h", type=float, default=4.0,
        help="tanh half-width of the transition, in units of the mean "
        "smoothing length.",
    )
    parser.add_argument(
        "--bulk-temperature-K", type=float, default=1000.0,
        help="Bulk gas temperature (K); requires SPH:initial_temperature: 0 "
        "in params.yml so the per-particle value survives.",
    )
    parser.add_argument(
        "--density-ratio", type=float, default=2.0,
        help="contrast variant only: density ratio across the shear layer, "
        "via a particle-mass ratio at uniform number density.",
    )
    parser.add_argument(
        "--source-geometry", type=str, default="blobs", choices=["blobs", "slab"],
        help="blobs: two Gaussian FUV/LW blobs, one per stream (gated). "
        "slab: a y-dependent slab seeded on the lower interface (report only).",
    )
    parser.add_argument(
        "--pulse-amplitude", type=float, default=1.0,
        help="Peak seeded FUVSpecificEnergy/LWSpecificEnergy (internal units).",
    )
    parser.add_argument(
        "--pulse-sigma-h", type=float, default=2.0,
        help="Gaussian width, in units of the mean smoothing length.",
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


########################################
# main
########################################

opt = parse_options()

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

# Mean density
rho = opt.rho  # atom/cc
rho = rho * constants.m_p / units.cm**3

# Gas particle mass
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
        N = int(np.ceil(M.to(units.Msun) / m))
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
print("Total mas                             : {}".format(M.to(units.Msun)))
print("Number of particles                   : {}".format(N))

# Convert to code units
m = m.to(UnitMass).value
L = L.to(UnitLength).value
rho = rho.to(UnitMass / UnitLength**3).value

# L_code is the box length, already converted to code units above.
L_code = L

# Pre-relaxed glass (getGlass.sh): random placement's shot noise would
# swamp the layer's own smooth velocity/density profile.
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

ids = np.arange(N)
h_mean = 1.2348 * L_code / N ** (1.0 / 3.0)  # 1.2348 = SPH:resolution_eta

print("Inter-particle distance (code unit)   : {}".format(L / N ** (1 / 3.0)))

#####################
# Velocity profile: KelvinHelmholtzGrowthRate_3D's two-interface profile,
# seed perturbation removed, tanh smoothing instead of a piecewise ramp.
#####################
y = pos[:, 1] / L_code  # in [0,1)
d = opt.layer_width_h * h_mean / L_code  # tanh half-width, box fraction
V = opt.v_shear_km_s / (UnitVelocity_in_cgs / 1e5)  # km/s -> code velocity
vx = 0.5 * V * (np.tanh((y - 0.25) / d) - np.tanh((y - 0.75) / d) - 1.0)
vel = np.zeros([N, 3])
vel[:, 0] = vx

print(f"variant                               : {opt.variant}")
print(f"v_shear (km/s)                        : {opt.v_shear_km_s}")
print(f"d/L                                   : {d}")
print(f"d/h                                   : {d * L_code / h_mean}")

#####################
# Masses and temperatures.
#####################
if opt.variant == "shear":
    mass = np.ones(N) * m
    T = np.full(N, opt.bulk_temperature_K)
    print(f"stream T (K)                          : {opt.bulk_temperature_K} (uniform)")
    print(f"stream mass (code)                    : {m} (uniform)")
else:
    R = opt.density_ratio
    s = 0.5 * (np.tanh((y - 0.25) / d) - np.tanh((y - 0.75) / d) + 1.0)  # 0 -> 1
    mass = m * (1.0 + (R - 1.0) * s)
    T = opt.bulk_temperature_K * (1.0 + (R - 1.0) * s) ** -1.0
    print(f"stream T (K)                          : {opt.bulk_temperature_K} / "
          f"{opt.bulk_temperature_K / R} (bulk / +V/2 stream)")
    print(f"stream mass (code)                    : {m} / {m * R} (bulk / +V/2 stream)")

UnitVelocity_cgs2 = UnitVelocity_in_cgs**2
u = np.array(
    [internal_energy_from_temperature_cgs(t) / UnitVelocity_cgs2 for t in T]
)

#####################
# Seeded FUV/LW field: two Gaussian blobs (gated) or an x-independent slab
# (report only). No star: see README.
#####################
u_fuv = np.zeros(N)
u_lw = np.zeros(N)
sigma = opt.pulse_sigma_h * h_mean
A = opt.pulse_amplitude

if opt.source_geometry == "blobs":
    blob_A_centre = np.array([0.25 * L_code, 0.50 * L_code, 0.50 * L_code])
    blob_B_centre = np.array([0.75 * L_code, 0.00 * L_code, 0.50 * L_code])

    def gaussian_blob(centre):
        dx = pos - centre
        dx -= L_code * np.round(dx / L_code)
        r2 = np.sum(dx**2, axis=1)
        return A * np.exp(-0.5 * r2 / sigma**2), dx

    u_A, dx_A = gaussian_blob(blob_A_centre)
    u_B, dx_B = gaussian_blob(blob_B_centre)
    u_fuv = u_A + u_B
    u_lw = u_fuv.copy()

    print(f"Blob A centre (code)                  : {blob_A_centre}")
    print(f"Blob B centre (code)                  : {blob_B_centre}")
    print(f"sigma/h                               : {sigma / h_mean}")

    r_A = np.sqrt(np.sum(dx_A**2, axis=1))
    r_B = np.sqrt(np.sum(dx_B**2, axis=1))
    in_A = r_A <= 3.0 * sigma
    in_B = r_B <= 3.0 * sigma
    N_A, N_B = int(in_A.sum()), int(in_B.sum())
    m_A, m_B = float(mass[in_A].sum()), float(mass[in_B].sum())
    mu_A = float((mass * u_fuv)[in_A].sum())
    mu_B = float((mass * u_fuv)[in_B].sum())
    imbalance_N = (N_A - N_B) / (N_A + N_B) if (N_A + N_B) > 0 else float("nan")
    imbalance_m = (m_A - m_B) / (m_A + m_B) if (m_A + m_B) > 0 else float("nan")
    imbalance_mu = (mu_A - mu_B) / (mu_A + mu_B) if (mu_A + mu_B) > 0 else float("nan")
    print(f"Blob A: N={N_A}  sum(m)={m_A:.6e}  sum(m*u)={mu_A:.6e}")
    print(f"Blob B: N={N_B}  sum(m)={m_B:.6e}  sum(m*u)={mu_B:.6e}")
    print(
        f"Population imbalance: (N_A-N_B)/(N_A+N_B)={imbalance_N:.4%}, "
        f"mass-weighted {imbalance_m:.4%}, (m*u)-weighted {imbalance_mu:.4%}"
    )
else:
    # slab: x-independent, seeded on the lower (y=0.25) interface.
    dy = (pos[:, 1] - 0.25 * L_code)
    dy -= L_code * np.round(dy / L_code)
    u_fuv = A * np.exp(-0.5 * dy**2 / sigma**2)
    u_lw = u_fuv.copy()
    print(f"Slab centred at y=0.25*L, sigma/h      : {sigma / h_mean}")

rho_arr = np.ones(N) * rho

#####################
# Finally write the ICs in the file
#####################

fileOutput = h5py.File(opt.outputfilename, "w")
print("{} saved.".format(opt.outputfilename))

grp = fileOutput.create_group("/Header")
grp.attrs["BoxSize"] = [L, L, L]
grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 0, 0]
grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
grp.attrs["Time"] = 0.0
grp.attrs["NumFileOutputsPerSnapshot"] = 1
grp.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
grp.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
grp.attrs["Dimension"] = 3

grp = fileOutput.create_group("/Units")
grp.attrs["Unit length in cgs (U_L)"] = UnitLength_in_cgs
grp.attrs["Unit mass in cgs (U_M)"] = UnitMass_in_cgs
grp.attrs["Unit time in cgs (U_t)"] = UnitTime_in_cgs
grp.attrs["Unit current in cgs (U_I)"] = UnitCurrent_in_cgs
grp.attrs["Unit temperature in cgs (U_T)"] = UnitTemp_in_cgs

grp = fileOutput.create_group("/PartType0")
grp.create_dataset("Coordinates", data=pos, dtype="d")
grp.create_dataset("Velocities", data=vel, dtype="f")
grp.create_dataset("Masses", data=mass, dtype="f")
grp.create_dataset("SmoothingLength", data=h, dtype="f")
grp.create_dataset("InternalEnergy", data=u, dtype="f")
grp.create_dataset("ParticleIDs", data=ids, dtype="L")
grp.create_dataset("Densities", data=rho_arr, dtype="f")
grp.create_dataset("FUVSpecificEnergy", data=u_fuv, dtype="f")
grp.create_dataset("LWSpecificEnergy", data=u_lw, dtype="f")

fileOutput.close()
