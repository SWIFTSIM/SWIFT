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
"""
Write a periodic glass box of two phases in pressure balance: x < L/2 cold
(temperature T, density n, particle mass m), x >= L/2 hot (16 T, n/16,
m/16). Both phases share the glass positions, hence the smoothing length;
the sound speed differs by 4, so the CFL time steps sit two bins apart.
One star sits in the hot phase near the seam at x = L/2.
"""

import argparse
import os

import h5py
import numpy as np
from astropy import constants, units

UNIT_MASS_CGS = 1.988409870698051e43
UNIT_LENGTH_CGS = 3.0856775814913673e21
UNIT_VELOCITY_CGS = 1e5
UNIT_TIME_CGS = UNIT_LENGTH_CGS / UNIT_VELOCITY_CGS
CONTRAST = 16.0


def parse_options() -> argparse.Namespace:
    """Parse the command line."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--level", type=int, default=5, help="N = (2**level)**3")
    parser.add_argument(
        "--rho", type=float, default=1e3, help="Cold density, atom/cm^3"
    )
    parser.add_argument("--mass", type=float, default=0.1, help="Cold gas mass, Msun")
    parser.add_argument("--T_cold", type=float, default=500.0, help="Cold phase, K")
    parser.add_argument("--star_mass", type=float, default=29.7, help="Msun")
    parser.add_argument(
        "--star_offset",
        type=float,
        default=2.0,
        help="Star distance from the seam at x = L/2, in interparticle "
        "spacings, on the hot side",
    )
    parser.add_argument("-o", dest="output", default="ICs_isrf_hierarchy.hdf5")
    return parser.parse_args()


def internal_energy(T: float, hydrogen_mass_fraction: float = 0.752) -> float:
    """Return the specific internal energy of neutral gas, internal units.

    Parameters
    ----------
    T : float
        Temperature in K, below the 1e4 K ionisation switch of SWIFT.
    hydrogen_mass_fraction : float
        Hydrogen mass fraction.

    Returns
    -------
    float
        u = k T / ((gamma - 1) mu m_p), gamma = 5/3.
    """
    if T >= 1e4:
        raise ValueError("Both phases must stay below 1e4 K (constant mu).")
    mu = 4.0 / (1.0 + 3.0 * hydrogen_mass_fraction)
    u_cgs = 1.380649e-16 * T * 1.5 / (mu * 1.67262192369e-24)
    return u_cgs / UNIT_VELOCITY_CGS**2


def main() -> None:
    """Write the initial conditions."""
    opt = parse_options()
    n_side = 2**opt.level
    n_gas = n_side**3

    rho = (opt.rho * constants.m_p / units.cm**3).to(units.g / units.cm**3).value
    m_cold = (opt.mass * units.Msun).to(units.g).value
    boxsize = (n_gas * m_cold / rho) ** (1.0 / 3.0) / UNIT_LENGTH_CGS
    m_cold /= UNIT_MASS_CGS

    name = f"glassCube_{n_side}.hdf5"
    if not os.path.exists(name):
        raise RuntimeError(f"Missing {name}. Run ./getGlass.sh {n_side}.")
    with h5py.File(name, "r") as glass:
        pos = np.mod(glass["/PartType0/Coordinates"][:, :], 1.0) * boxsize
        h = glass["/PartType0/SmoothingLength"][:] * boxsize
    if pos.shape[0] != n_gas:
        raise RuntimeError(f"{name} has {pos.shape[0]} particles, not {n_gas}.")

    hot = pos[:, 0] >= 0.5 * boxsize
    mass = np.where(hot, m_cold / CONTRAST, m_cold)
    u = np.where(
        hot,
        internal_energy(CONTRAST * opt.T_cold),
        internal_energy(opt.T_cold),
    )

    spacing = boxsize / n_side
    pos_star = np.array(
        [[0.5 * boxsize + opt.star_offset * spacing, 0.5 * boxsize, 0.5 * boxsize]]
    )
    m_star = (opt.star_mass * units.Msun).to(units.g).value / UNIT_MASS_CGS

    print(
        f"Box size: {boxsize:.6e} kpc, cold gas: {int((~hot).sum())}, "
        f"hot gas: {int(hot.sum())}, T = {opt.T_cold} K and {CONTRAST * opt.T_cold} K"
    )

    with h5py.File(opt.output, "w") as f:
        grp = f.create_group("/Header")
        grp.attrs["BoxSize"] = [boxsize] * 3
        grp.attrs["NumPart_Total"] = [n_gas, 0, 0, 0, 1, 0]
        grp.attrs["NumPart_Total_HighWord"] = [0] * 6
        grp.attrs["NumPart_ThisFile"] = [n_gas, 0, 0, 0, 1, 0]
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
        grp.create_dataset("Masses", data=mass, dtype="f")
        grp.create_dataset("SmoothingLength", data=h, dtype="f")
        grp.create_dataset("InternalEnergy", data=u, dtype="f")
        grp.create_dataset("ParticleIDs", data=np.arange(1, n_gas + 1), dtype="L")

        grp = f.create_group("/PartType4")
        grp.create_dataset("Coordinates", data=pos_star, dtype="d")
        grp.create_dataset("Velocities", data=np.zeros((1, 3)), dtype="f")
        grp.create_dataset("Masses", data=[m_star], dtype="f")
        grp.create_dataset("BirthMass", data=[m_star], dtype="f")
        grp.create_dataset("BirthTime", data=[0.0], dtype="f")
        grp.create_dataset("ParticleIDs", data=[n_gas + 1], dtype="L")
        grp.create_dataset("SmoothingLength", data=[3.0 * spacing], dtype="f")
        grp.create_dataset("StellarParticleType", data=[0], dtype="i")


if __name__ == "__main__":
    main()
