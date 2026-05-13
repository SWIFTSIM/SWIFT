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
args = parser.parse_args()
sidm_fname = args.sidm

lw = 1.2

S_IN_GYR = 3.1536e16
G = 4.302e-6  # kpc (km/s)^2 / Msun
rho_crit = 127.0  # Msun / kpc^3  (z=0)
KM_S_TO_KPC_PER_GYR = 1.022
MSUN_IN_G = 1.988409e33
KPC_IN_CM = 3.085677e21

# ---- sigma/m: 100 cm^2/g -> kpc^2/Msun ----
sigma_over_m = 100.0 / KPC_IN_CM**2 * MSUN_IN_G


# ------------------ Units ------------------
def get_unit_time(f):
    return f["Units"].attrs["Unit time in cgs (U_t)"][0]


def density_conversion(f):
    u_mass = f["Units"].attrs["Unit mass in cgs (U_M)"][0]
    u_length = f["Units"].attrs["Unit length in cgs (U_L)"][0]
    return (u_mass / MSUN_IN_G) / (u_length / KPC_IN_CM) ** 3


# ------------------ Helpers ------------------
def compute_radius(pos):
    center = np.median(pos, axis=0)
    return np.sqrt(((pos - center) ** 2).sum(axis=1))


# ------------------ NFW ------------------
def nfw_params(M200, c):
    r200 = (3 * M200 / (4 * np.pi * 200 * rho_crit)) ** (1 / 3)
    rs = r200 / c
    f_c = np.log(1 + c) - c / (1 + c)
    rho_s = M200 / (4 * np.pi * rs**3 * f_c)
    return r200, rs, rho_s


def nfw_density(r, rho_s, rs):
    x = r / rs
    return rho_s / (x * (1 + x) ** 2)


def nfw_mass_enclosed(r, rho_s, rs):
    x = r / rs
    return 4 * np.pi * rho_s * rs**3 * (np.log(1 + x) - x / (1 + x))


def scattering_rate_nfw(r, M200, c, sigma_over_m):
    _, rs, rho_s = nfw_params(M200, c)
    rho = nfw_density(r, rho_s, rs)
    Menc = nfw_mass_enclosed(r, rho_s, rs)
    sigma_1d = np.sqrt(G * Menc / np.maximum(r, 1e-10)) / np.sqrt(2)
    v_rel = 4 * sigma_1d / np.sqrt(np.pi)  # km/s
    return rho * sigma_over_m * v_rel * KM_S_TO_KPC_PER_GYR


# ------------------ Radial binning ------------------
def median_rate_vs_r(df, cols, bins):
    df = df.copy()
    df["r_bin"] = pd.cut(df["r"], bins=bins)
    grouped = df.groupby("r_bin", observed=True)[cols].median()
    r_mid = 0.5 * (bins[:-1] + bins[1:])
    # align r_mid to the bins that groupby kept
    bin_labels = grouped.index
    r_mid_out = np.array([0.5 * (b.left + b.right) for b in bin_labels])
    return r_mid_out, grouped


def main():

    with h5py.File(sidm_fname, "r") as f:
        pos = np.array(f["SIDMParticles"]["Coordinates"])
        rates = np.array(f["SIDMParticles"]["Rates"])
        vel = np.array(f["SIDMParticles"]["Velocities"])
        rho = np.array(f["SIDMParticles"]["Densities"]) * density_conversion(f)
        u_time = get_unit_time(f)

    ## Calculate scattering rate based on particle densities and velocities
    v_rel_kpc_gyr = np.sqrt((vel**2).sum(axis=1)) * KM_S_TO_KPC_PER_GYR
    rate_particle = rho * sigma_over_m * v_rel_kpc_gyr

    df_sidm = pd.DataFrame(
        {
            "rate_sim": rates / u_time * S_IN_GYR,
            "rate_part": rate_particle,
            "r": compute_radius(pos),
        }
    )

    # ---- Radial bins ----
    bins = np.logspace(-2, 1.7, 40)
    r_mid, grp = median_rate_vs_r(df_sidm, ["rate_sim", "rate_part"], bins)

    sim_med = grp["rate_sim"].values
    part_med = grp["rate_part"].values

    # ---- Analytic NFW ----
    r200, rs, rho_s = nfw_params(M200=5e9, c=12)
    analytic = scattering_rate_nfw(r_mid, M200=5e9, c=12, sigma_over_m=sigma_over_m)

    # ---- Plot ----
    fig, ax = plt.subplots(1, 1)
    ax.plot(r_mid, sim_med, lw=lw, label="Stored rates")
    ax.plot(
        r_mid, part_med, lw=lw, label=r"$\rho\,(\sigma/m)\,|v|$ (particles)", ls="--"
    )
    ax.plot(r_mid, analytic, lw=lw, label="Analytic NFW")
    ax.axvline(rs, color="gray", ls=":", lw=0.8, label=f"$r_s$={rs:.2f} kpc")
    ax.axvline(r200, color="gray", ls="--", lw=0.8, label=f"$r_{{200}}$={r200:.1f} kpc")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("r [kpc]")
    ax.set_ylabel(r"Scattering rate [Gyr$^{-1}$]")
    ax.legend()
    fig.tight_layout()
    plt.savefig("scattering_rates.png", dpi=200)


if __name__ == "__main__":
    main()
