###############################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Katy Proctor (katy.proctor@fysik.su.se)
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
#################################################################################

# Uniform-density periodic box of SIDM particles for testing the scattering
# rate in a cosmological run.
#
# Positions are uniform random (Poisson) in a comoving box of side L, and the
# particle mass is set so that the comoving density is Omega_cdm rho_crit,0
# (the run uses Omega_b = 0). Peculiar velocities are Maxwellian with 1D
# dispersion sigma_1D and zero mean. SWIFT reads the IC velocities as peculiar
# velocities. The expected rate is then
#
#   Gamma = (rho_comoving / a^3) (sigma/m) 4 sigma_1D / sqrt(pi).

import argparse
import h5py
import numpy as np

# Units of the ICs (and of params.yml)
UNIT_MASS_CGS = 1.98841e43  # 10^10 Msun
UNIT_LENGTH_CGS = 3.08567758e24  # Mpc
UNIT_VELOCITY_CGS = 1e5  # km/s

G_CGS = 6.67430e-8
MPC_CGS = 3.08567758e24

parser = argparse.ArgumentParser()
parser.add_argument("-n", type=int, default=100, help="Particles per dimension (N = n^3)")
parser.add_argument("-L", type=float, default=10.0, help="Comoving box size [Mpc]")
parser.add_argument("--sigma", type=float, default=100.0, help="1D velocity dispersion [km/s]")
parser.add_argument("--h-hubble", type=float, default=0.6777, help="Reduced Hubble constant")
parser.add_argument("--omega-cdm", type=float, default=0.3, help="Omega_cdm")
parser.add_argument("--eta", type=float, default=1.2348, help="Initial guess for h")
parser.add_argument("--seed", type=int, default=42)
parser.add_argument("-o", "--output", default="uniform_box.hdf5")


def main():
    args = parser.parse_args()
    rng = np.random.default_rng(args.seed)
    n_part = args.n**3

    # Comoving density in internal units
    H0 = 100 * args.h_hubble * 1e5 / MPC_CGS  # 1/s
    rho_crit_cgs = 3 * H0**2 / (8 * np.pi * G_CGS)
    rho_crit = rho_crit_cgs / UNIT_MASS_CGS * UNIT_LENGTH_CGS**3
    mass = args.omega_cdm * rho_crit * args.L**3 / n_part

    pos = rng.random((n_part, 3)) * args.L
    vel = rng.normal(0.0, args.sigma, (n_part, 3))
    vel -= vel.mean(axis=0)
    h = args.eta * args.L / args.n

    with h5py.File(args.output, "w") as f:
        num = np.zeros(8, dtype=np.uint32)
        num[7] = n_part
        hdr = f.create_group("Header")
        hdr.attrs["BoxSize"] = [args.L] * 3
        hdr.attrs["NumPart_ThisFile"] = num
        hdr.attrs["NumPart_Total"] = num
        hdr.attrs["NumPart_Total_HighWord"] = np.zeros(8, dtype=np.uint32)
        hdr.attrs["MassTable"] = np.zeros(8)
        hdr.attrs["Time"] = 0.0
        hdr.attrs["NumFilesPerSnapshot"] = 1
        hdr.attrs["Flag_Entropy_ICs"] = 0
        hdr.attrs["Dimension"] = 3

        units = f.create_group("Units")
        units.attrs["Unit mass in cgs (U_M)"] = UNIT_MASS_CGS
        units.attrs["Unit length in cgs (U_L)"] = UNIT_LENGTH_CGS
        units.attrs["Unit time in cgs (U_t)"] = UNIT_LENGTH_CGS / UNIT_VELOCITY_CGS
        units.attrs["Unit current in cgs (U_I)"] = 1.0
        units.attrs["Unit temperature in cgs (U_T)"] = 1.0

        g = f.create_group("PartType7")
        g.create_dataset("Coordinates", data=pos, dtype="d")
        g.create_dataset("Velocities", data=vel, dtype="f")
        g.create_dataset("Masses", data=np.full(n_part, mass), dtype="f")
        g.create_dataset("SmoothingLength", data=np.full(n_part, h), dtype="f")
        g.create_dataset("ParticleIDs", data=np.arange(1, n_part + 1), dtype="L")

    print(f"Wrote {n_part} particles to {args.output}")
    print(f"  L = {args.L} Mpc, m = {mass * 1e10:.4e} Msun, "
          f"rho_comoving = {args.omega_cdm * rho_crit * 1e10 / 1e9:.4e} Msun/kpc^3")
    print(f"  sigma_1D = {args.sigma} km/s")


if __name__ == "__main__":
    main()
