################################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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
"""Write a uniform glass box with an optional star and a seeded ISRF.

Every input is PHYSICAL at the starting scale factor ``a = 1/(1+z)``. SWIFT
reads cosmological positions and densities as comoving, and converts the IC
internal energy from physical to comoving itself, so the file holds

    comoving box size          = L / a
    comoving density           = rho * a**3
    internal energy            = u   (physical)
    star BirthTime             = a   (read as the birth scale factor)

SWIFT stores the birth scale factor in single precision. It is rounded down to
the nearest float32 value, otherwise the star is born up to one float32 step in
a after the start (7.5e-9 in a at a = 0.1, a quarter of the H2 run's span).
The seeded ``PESpecificEnergy`` and ``LWSpecificEnergy`` are physical and
mass-specific, so they carry no scale factor. With ``--redshift 0`` the file is
a plain non-cosmological IC.
"""

import argparse
import os

import h5py
import numpy as np
import yaml
from astropy import constants
from astropy import units
from cosmo_timeline import hubble_rate

UNIT_MASS_CGS = 1.988409870698051e43
UNIT_LENGTH_CGS = 3.0856775814913673e21
UNIT_VELOCITY_CGS = 1e5
UNIT_TIME_CGS = UNIT_LENGTH_CGS / UNIT_VELOCITY_CGS
GAMMA = 5.0 / 3.0
HYDROGEN_MASS_FRACTION = 0.76


def parse_options() -> argparse.Namespace:
    """Parse the command line."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level", type=int, default=5, help="N = (2**level)**3")
    parser.add_argument(
        "--rho", type=float, default=1.0, help="Physical density, atom/cm^3"
    )
    parser.add_argument(
        "--mass", type=float, default=1.0, help="Gas particle mass, Msun"
    )
    parser.add_argument(
        "--temperature", type=float, default=100.0, help="Physical temperature, K"
    )
    parser.add_argument(
        "--redshift",
        type=float,
        default=0.0,
        help="Starting redshift, 0 = no cosmology",
    )
    parser.add_argument(
        "--u-pe",
        type=float,
        default=0.0,
        help="Seeded physical PE specific energy, erg/g",
    )
    parser.add_argument(
        "--u-lw",
        type=float,
        default=0.0,
        help="Seeded physical LW specific energy, erg/g",
    )
    parser.add_argument(
        "--star-mass",
        type=float,
        default=0.0,
        help="Mass of the central star in Msun, 0 = no star",
    )
    parser.add_argument(
        "--star-age",
        type=float,
        default=0.0,
        help="Age of the star at the start, internal time units. Without "
        "cosmology the run must start at this time (run.sh does so)",
    )
    parser.add_argument(
        "--params", default="params.yml", help="Parameter file, for the cosmology"
    )
    parser.add_argument("-o", dest="output", default="ICs_isrf_cosmology.hdf5")
    return parser.parse_args()


def main() -> None:
    """Write the IC file."""
    opt = parse_options()
    a = 1.0 / (1.0 + opt.redshift)

    n_side = 2**opt.level
    n_gas = n_side**3
    mass_msun = opt.mass
    rho_cgs = opt.rho * constants.m_p.cgs.value
    total_mass_cgs = n_gas * mass_msun * constants.M_sun.cgs.value
    box_physical_cgs = (total_mass_cgs / rho_cgs) ** (1.0 / 3.0)

    # Neutral atomic gas, the mean molecular weight SWIFT also uses below
    # its ionization temperature.
    mu = 4.0 / (1.0 + 3.0 * HYDROGEN_MASS_FRACTION)
    u_physical_cgs = (
        constants.k_B.cgs.value
        * opt.temperature
        / ((GAMMA - 1.0) * mu * constants.m_p.cgs.value)
    )

    box = box_physical_cgs / a / UNIT_LENGTH_CGS
    rho_comoving = rho_cgs * a**3 * UNIT_LENGTH_CGS**3 / UNIT_MASS_CGS
    u_internal = u_physical_cgs / UNIT_VELOCITY_CGS**2
    mass = mass_msun * constants.M_sun.cgs.value / UNIT_MASS_CGS
    energy_unit = UNIT_VELOCITY_CGS**2

    glass_filename = f"glassCube_{n_side}.hdf5"
    if not os.path.exists(glass_filename):
        raise RuntimeError(f"Missing {glass_filename}. Run ./getGlass.sh {n_side}.")
    with h5py.File(glass_filename, "r") as glass:
        pos = glass["/PartType0/Coordinates"][:, :]
        h_glass = glass["/PartType0/SmoothingLength"][:]
    if pos.shape[0] != n_gas:
        raise RuntimeError(f"{glass_filename} has {pos.shape[0]} particles.")
    pos = (pos - pos.min()) / (pos.max() - pos.min() + 1e-6) * box
    h = h_glass * box

    print(f"Scale factor               : {a:.8g}")
    print(f"Physical box size          : {box * a * 1e3:.6g} pc")
    print(f"Comoving box size          : {box:.6g} kpc")
    recovered_rho = rho_comoving / a**3 * UNIT_MASS_CGS / UNIT_LENGTH_CGS**3
    print(
        "Physical density (check)   : "
        f"{recovered_rho / constants.m_p.cgs.value:.6g} atom/cm^3"
    )
    recovered_T = (
        u_internal
        * UNIT_VELOCITY_CGS**2
        * (GAMMA - 1.0)
        * mu
        * constants.m_p.cgs.value
        / constants.k_B.cgs.value
    )
    print(f"Physical temperature (check): {recovered_T:.6g} K")
    print(f"Gas particles              : {n_gas}")

    with h5py.File(opt.output, "w") as handle:
        n_star = 1 if opt.star_mass > 0.0 else 0
        header = handle.create_group("/Header")
        header.attrs["BoxSize"] = [box, box, box]
        header.attrs["NumPart_Total"] = [n_gas, 0, 0, 0, n_star, 0]
        header.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
        header.attrs["NumPart_ThisFile"] = [n_gas, 0, 0, 0, n_star, 0]
        header.attrs["Time"] = 0.0
        header.attrs["NumFileOutputsPerSnapshot"] = 1
        header.attrs["MassTable"] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        header.attrs["Flag_Entropy_ICs"] = [0, 0, 0, 0, 0, 0]
        header.attrs["Dimension"] = 3

        unit_group = handle.create_group("/Units")
        unit_group.attrs["Unit length in cgs (U_L)"] = UNIT_LENGTH_CGS
        unit_group.attrs["Unit mass in cgs (U_M)"] = UNIT_MASS_CGS
        unit_group.attrs["Unit time in cgs (U_t)"] = UNIT_TIME_CGS
        unit_group.attrs["Unit current in cgs (U_I)"] = 1.0
        unit_group.attrs["Unit temperature in cgs (U_T)"] = 1.0

        gas = handle.create_group("/PartType0")
        gas.create_dataset("Coordinates", data=pos, dtype="d")
        gas.create_dataset("Velocities", data=np.zeros((n_gas, 3)), dtype="f")
        gas.create_dataset("Masses", data=np.full(n_gas, mass), dtype="f")
        gas.create_dataset("SmoothingLength", data=h, dtype="f")
        gas.create_dataset("InternalEnergy", data=np.full(n_gas, u_internal), dtype="f")
        gas.create_dataset("ParticleIDs", data=np.arange(n_gas), dtype="L")
        gas.create_dataset("Density", data=np.full(n_gas, rho_comoving), dtype="f")
        gas.create_dataset(
            "PESpecificEnergy",
            data=np.full(n_gas, opt.u_pe / energy_unit),
            dtype="f",
        )
        gas.create_dataset(
            "LWSpecificEnergy",
            data=np.full(n_gas, opt.u_lw / energy_unit),
            dtype="f",
        )

        if n_star:
            star_mass = opt.star_mass * constants.M_sun.cgs.value / UNIT_MASS_CGS
            star = handle.create_group("/PartType4")
            star.create_dataset("Coordinates", data=np.full((1, 3), box / 2), dtype="d")
            star.create_dataset("Velocities", data=np.zeros((1, 3)), dtype="f")
            star.create_dataset("Masses", data=[star_mass], dtype="f")
            star.create_dataset("ParticleIDs", data=[n_gas + 1], dtype="L")
            star.create_dataset("SmoothingLength", data=[3 * box / n_side], dtype="f")
            star.create_dataset("BirthMass", data=[star_mass], dtype="f")
            # GEAR skips feedback for a negative birth time; without cosmology
            # run.sh starts the run at time star_age instead.
            birth = 0.0
            if opt.redshift > 0.0:
                with open(opt.params) as handle:
                    cosmology = yaml.safe_load(handle)["Cosmology"]
                # H t << 1 over the age, so ln a is linear in t.
                a_birth = a * np.exp(-hubble_rate(a, cosmology) * opt.star_age)
                birth = np.float32(a_birth)
                if float(birth) > a_birth:
                    birth = np.nextafter(birth, np.float32(0.0))
                print(
                    f"Star birth scale factor    : {float(birth):.10g} "
                    f"(start {a:.10g})"
                )
            star.create_dataset("BirthTime", data=[birth], dtype="f")
            star.create_dataset("StellarParticleType", data=[0], dtype="i")

    print(f"{opt.output} saved.")


if __name__ == "__main__":
    main()
