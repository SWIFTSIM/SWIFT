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
"""Write a uniform periodic glass box with one, two or a lattice of stars."""

import argparse
import os

import h5py
import numpy as np
from astropy import constants, units

UNIT_MASS_CGS = 1.988409870698051e43
UNIT_LENGTH_CGS = 3.0856775814913673e21
UNIT_VELOCITY_CGS = 1e5
UNIT_TIME_CGS = UNIT_LENGTH_CGS / UNIT_VELOCITY_CGS


def parse_options() -> argparse.Namespace:
    """Parse the command line."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level", type=int, default=5, help="N = (2**level)**3")
    parser.add_argument("--rho", type=float, default=1e3, help="Density, atom/cm^3")
    parser.add_argument("--mass", type=float, default=0.1, help="Gas mass, Msun")
    parser.add_argument("--star_mass", type=float, default=29.7, help="Msun")
    parser.add_argument(
        "--sources",
        choices=["A", "B", "AB", "lattice", "lattice_single"],
        default="AB",
        help="A: one star at (1/2 - s/2, 1/2, 1/2) L; B: one star at (1/2 + "
        "s/2, 1/2, 1/2) L, s the separation; AB: both; lattice: n_side^3 stars at the cell centres of a "
        "cubic lattice of spacing L/n_side; lattice_single: the first star "
        "of that lattice alone.",
    )
    parser.add_argument("--n_side", type=int, default=8, help="Lattice stars per side")
    parser.add_argument(
        "--separation", type=float, default=0.25, help="A to B distance over L"
    )
    parser.add_argument("-o", dest="output", default="ICs_isrf_multi_source.hdf5")
    return parser.parse_args()


def star_positions(
    sources: str, n_side: int, separation: float, boxsize: float
) -> np.ndarray:
    """Return the star positions for a source configuration.

    Parameters
    ----------
    sources : str
        Configuration name (see ``--sources``).
    n_side : int
        Lattice stars per side.
    separation : float
        Distance between stars A and B over the box size.
    boxsize : float
        Box side length, internal units.

    Returns
    -------
    np.ndarray
        Positions, shape (n_stars, 3).
    """
    a = np.array([0.5 - 0.5 * separation, 0.5, 0.5])
    b = np.array([0.5 + 0.5 * separation, 0.5, 0.5])
    if sources == "A":
        frac = a[None, :]
    elif sources == "B":
        frac = b[None, :]
    elif sources == "AB":
        frac = np.vstack([a, b])
    else:
        c = (np.arange(n_side) + 0.5) / n_side
        frac = np.stack(np.meshgrid(c, c, c, indexing="ij"), axis=-1).reshape(-1, 3)
        if sources == "lattice_single":
            frac = frac[:1]
    return frac * boxsize


def main() -> None:
    """Write the initial conditions."""
    opt = parse_options()
    n_side_gas = 2**opt.level
    n_gas = n_side_gas**3

    rho = (opt.rho * constants.m_p / units.cm**3).to(units.g / units.cm**3).value
    m_gas = (opt.mass * units.Msun).to(units.g).value
    boxsize = (n_gas * m_gas / rho) ** (1.0 / 3.0) / UNIT_LENGTH_CGS
    m_gas /= UNIT_MASS_CGS

    glass_file = f"glassCube_{n_side_gas}.hdf5"
    if not os.path.exists(glass_file):
        raise RuntimeError(f"Missing {glass_file}. Run ./getGlass.sh {n_side_gas}.")
    with h5py.File(glass_file, "r") as glass:
        pos = glass["/PartType0/Coordinates"][:, :] * boxsize
        h = glass["/PartType0/SmoothingLength"][:] * boxsize
    if pos.shape[0] != n_gas:
        raise RuntimeError(f"{glass_file} has {pos.shape[0]} particles, not {n_gas}.")
    pos = np.mod(pos, boxsize)

    pos_star = star_positions(opt.sources, opt.n_side, opt.separation, boxsize)
    n_star = pos_star.shape[0]
    m_star = (opt.star_mass * units.Msun).to(units.g).value / UNIT_MASS_CGS

    print(f"Box size: {boxsize:.6e} kpc, gas particles: {n_gas}, stars: {n_star}")

    with h5py.File(opt.output, "w") as f:
        grp = f.create_group("/Header")
        grp.attrs["BoxSize"] = [boxsize] * 3
        grp.attrs["NumPart_Total"] = [n_gas, 0, 0, 0, n_star, 0]
        grp.attrs["NumPart_Total_HighWord"] = [0] * 6
        grp.attrs["NumPart_ThisFile"] = [n_gas, 0, 0, 0, n_star, 0]
        grp.attrs["Time"] = 0.0
        grp.attrs["NumFileOutputsPerSnapshot"] = 1
        grp.attrs["MassTable"] = [0.0] * 6
        grp.attrs["Flag_Entropy_ICs"] = [0] * 6
        grp.attrs["Dimension"] = 3

        grp = f.create_group("/Units")
        grp.attrs["Unit length in cgs (U_L)"] = UNIT_LENGTH_CGS
        grp.attrs["Unit mass in cgs (U_M)"] = UNIT_MASS_CGS
        grp.attrs["Unit time in cgs (U_t)"] = UNIT_TIME_CGS
        grp.attrs["Unit current in cgs (U_I)"] = 1.0
        grp.attrs["Unit temperature in cgs (U_T)"] = 1.0

        grp = f.create_group("/PartType0")
        grp.create_dataset("Coordinates", data=pos, dtype="d")
        grp.create_dataset("Velocities", data=np.zeros((n_gas, 3)), dtype="f")
        grp.create_dataset("Masses", data=np.full(n_gas, m_gas), dtype="f")
        grp.create_dataset("SmoothingLength", data=h, dtype="f")
        grp.create_dataset("InternalEnergy", data=np.ones(n_gas), dtype="f")
        grp.create_dataset("ParticleIDs", data=np.arange(1, n_gas + 1), dtype="L")

        grp = f.create_group("/PartType4")
        grp.create_dataset("Coordinates", data=pos_star, dtype="d")
        grp.create_dataset("Velocities", data=np.zeros((n_star, 3)), dtype="f")
        grp.create_dataset("Masses", data=np.full(n_star, m_star), dtype="f")
        grp.create_dataset("BirthMass", data=np.full(n_star, m_star), dtype="f")
        grp.create_dataset("BirthTime", data=np.zeros(n_star), dtype="f")
        grp.create_dataset(
            "ParticleIDs", data=np.arange(n_gas + 1, n_gas + n_star + 1), dtype="L"
        )
        grp.create_dataset(
            "SmoothingLength",
            data=np.full(n_star, 3.0 * boxsize / n_side_gas),
            dtype="f",
        )
        grp.create_dataset("StellarParticleType", data=np.zeros(n_star), dtype="i")


if __name__ == "__main__":
    main()
