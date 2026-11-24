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
"""Star-free periodic gas box with a uniform far-UV field written into the ICs."""

import argparse
import os

import h5py
import numpy as np

UNIT_MASS_CGS = 1.988409870698051e43  # 1e10 Msun
UNIT_LENGTH_CGS = 3.0856775814913673e21  # kpc
UNIT_VELOCITY_CGS = 1e5  # km/s
UNIT_TIME_CGS = UNIT_LENGTH_CGS / UNIT_VELOCITY_CGS
MSUN_CGS = 1.988409870698051e33
PROTON_MASS_CGS = 1.67262192369e-24
HYDROGEN_MASS_FRACTION = 0.76
SPEED_OF_LIGHT_CGS = 2.99792458e10
HABING_FLUX_CGS = 1.6e-3  # erg s^-1 cm^-2, G0 = 1


def parse_options():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--rho", type=float, default=10.0, help="Gas density in H atoms/cm3"
    )
    parser.add_argument(
        "--mass", type=float, default=0.1, help="Gas particle mass in Msun"
    )
    parser.add_argument(
        "--level", type=int, default=5, help="Resolution level: N = (2**level)**3"
    )
    parser.add_argument(
        "--g0",
        type=float,
        default=1.65,
        help="Far-UV field in Habing units (1.65 stands for chi = 1 of the Draine field, approximately)",
    )
    parser.add_argument(
        "--pe_fraction",
        type=float,
        default=1.0,
        help="Share of the field in the photoelectric band (the rest is in the Lyman-Werner band)",
    )
    parser.add_argument(
        "-o", dest="outputfilename", default="box.hdf5", help="Output file"
    )
    return parser.parse_args()


opt = parse_options()

N = (2**opt.level) ** 3
rho_cgs = opt.rho * PROTON_MASS_CGS / HYDROGEN_MASS_FRACTION
m_cgs = opt.mass * MSUN_CGS
L_cgs = (N * m_cgs / rho_cgs) ** (1.0 / 3.0)

glass_file = f"glassCube_{2**opt.level}.hdf5"
if not os.path.exists(glass_file):
    raise RuntimeError(f"Missing {glass_file}. Run ./getGlass.sh {2**opt.level} first.")
with h5py.File(glass_file, "r") as glass:
    pos = glass["/PartType0/Coordinates"][:, :]
    h_glass = glass["/PartType0/SmoothingLength"][:]
if pos.shape[0] != N:
    raise RuntimeError(f"{glass_file} has {pos.shape[0]} particles, expected {N}.")

L = L_cgs / UNIT_LENGTH_CGS
pos = (pos - pos.min()) / (pos.max() - pos.min() + 1e-6) * L
h = h_glass * L
mass = np.ones(N) * m_cgs / UNIT_MASS_CGS
rho = np.ones(N) * rho_cgs / (UNIT_MASS_CGS / UNIT_LENGTH_CGS**3)

# G0 = c rho (u_PE + u_LW) / F_Habing, with u the field energy per unit gas mass.
# The field is written for the nominal density: a particle whose SPH density differs
# from it sees a G0 that follows its own density.
u_field_cgs = opt.g0 * HABING_FLUX_CGS / (SPEED_OF_LIGHT_CGS * rho_cgs)
u_field = u_field_cgs / UNIT_VELOCITY_CGS**2
u_pe = np.ones(N) * u_field * opt.pe_fraction
u_lw = np.ones(N) * u_field * (1.0 - opt.pe_fraction)

print(f"Boxsize [pc]                : {L_cgs / 3.0856775814913673e18:.5g}")
print(f"Number of particles         : {N}")
print(f"Field energy per mass [erg/g]: {u_field_cgs:.6e}")

with h5py.File(opt.outputfilename, "w") as f:
    grp = f.create_group("/Header")
    grp.attrs["BoxSize"] = [L, L, L]
    grp.attrs["NumPart_Total"] = [N, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_Total_HighWord"] = [0, 0, 0, 0, 0, 0]
    grp.attrs["NumPart_ThisFile"] = [N, 0, 0, 0, 0, 0]
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
    grp.create_dataset("Velocities", data=np.zeros((N, 3)), dtype="f")
    grp.create_dataset("Masses", data=mass, dtype="f")
    grp.create_dataset("SmoothingLength", data=h, dtype="f")
    grp.create_dataset("InternalEnergy", data=np.ones(N), dtype="f")
    grp.create_dataset("ParticleIDs", data=np.arange(N), dtype="L")
    grp.create_dataset("Densities", data=rho, dtype="f")
    grp.create_dataset("PESpecificEnergy", data=u_pe, dtype="d")
    grp.create_dataset("LWSpecificEnergy", data=u_lw, dtype="d")
