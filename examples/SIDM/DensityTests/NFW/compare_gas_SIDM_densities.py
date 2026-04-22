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
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import h5py

parser = argparse.ArgumentParser()
parser.add_argument(
    "--sidm", type=str, help="SIDM snapshot file", default="snap/snapshot_0000.hdf5"
)
parser.add_argument(
    "--gas", type=str, help="Gas snapshot file", default="snap_gas/snapshot_0000.hdf5"
)

args = parser.parse_args()
sidm_fname = args.sidm
gas_fname = args.gas

MSUN_IN_G = 1.988409e33
KPC_IN_CM = 3.085677e21


def get_unit_conversions(f):
    u_mass = f["Units"].attrs["Unit mass in cgs (U_M)"][0]
    u_length = f["Units"].attrs["Unit length in cgs (U_L)"][0]
    return (u_mass / MSUN_IN_G) / (u_length / KPC_IN_CM) ** 3


def main():

    with h5py.File(sidm_fname, "r") as f:
        sidm_pid = np.array(f["SIDMParticles"]["ParticleIDs"])
        sidm_rho_sph = np.array(f["SIDMParticles"]["Densities"]) * get_unit_conversions(
            f
        )
    df_sidm = pd.DataFrame({"ParticleIDs": sidm_pid, "Densities": sidm_rho_sph})

    with h5py.File(gas_fname, "r") as f:
        gas_pid = np.array(f["GasParticles"]["ParticleIDs"])
        gas_rho_sph = np.array(f["GasParticles"]["Densities"]) * get_unit_conversions(f)
    df_gas = pd.DataFrame({"ParticleIDs": gas_pid, "Densities": gas_rho_sph})

    # match on PID
    df = df_sidm.merge(df_gas, on="ParticleIDs", how="left", suffixes=["_sidm", "_gas"])
    print("correct match: ", df.shape[0] == df_sidm.shape[0] == df_gas.shape[0])

    df["rho_check"] = df["Densities_sidm"] / df["Densities_gas"]

    print(df["rho_check"].min())
    print(df["rho_check"].max())

    fig, ax = plt.subplots(1, 1, figsize=(6.6, 6))
    ax.scatter(
        np.log10(df["Densities_sidm"]),
        np.log10(df["rho_check"]),
        s=1,
        alpha=0.1,
        lw=0.1,
    )
    ax.set_ylabel(r"$\log\, \rho_{\rm{SIDM}}/\rho_{\rm{gas}}$")
    ax.set_xlabel(r"$\log\, \rho_{\rm{SIDM}}$")
    plt.tight_layout()
    plt.savefig("gas_SIDM_density_comparison.png", dpi=250)


if __name__ == "__main__":
    main()
