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

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as mpe
import argparse
import h5py

parser = argparse.ArgumentParser()
parser.add_argument("file", type=str, help="Snapshot file")
parser.add_argument("--notex", action="store_true", help="Disable LaTeX")
parser.add_argument("-Rmin", type=float, default=0.005)
parser.add_argument("-Rmax", type=float, default=50)
parser.add_argument("-nbsamples", type=int, default=64)
parser.add_argument("-shift", type=float, default=2.0)

args = parser.parse_args()
fname = args.file

plt.rcParams.update(
    {
        "text.usetex": not args.notex,
        "font.size": 14 if not args.notex else 12,
        "font.family": "serif" if not args.notex else None,
    }
)


MSUN_IN_G = 1.988409e33
KPC_IN_CM = 3.085677e21

rsp = np.logspace(np.log10(args.Rmin), np.log10(args.Rmax), args.nbsamples)


def get_unit_conversions(f):
    u_mass = f["Units"].attrs["Unit mass in cgs (U_M)"][0]
    u_length = f["Units"].attrs["Unit length in cgs (U_L)"][0]
    return (u_mass / MSUN_IN_G) / (u_length / KPC_IN_CM) ** 3


def compute_radius(pos):
    pos -= np.median(pos, axis=0)
    return np.sqrt((pos**2).sum(axis=1))


def bin_centers(bins):
    # Geometric mid-points, appropriate for logarithmically spaced bins
    return np.sqrt(bins[:-1] * bins[1:])


def binned_density(r, mass, bins):
    m_cum = np.array([(mass[r < R]).sum() for R in bins])
    dM = np.diff(m_cum)
    shell_vol = (4.0 / 3.0) * np.pi * (bins[1:] ** 3 - bins[:-1] ** 3)
    return bin_centers(bins), dM / shell_vol


def sph_density_profile(r, rho_sph, bins):
    centers = bin_centers(bins)
    rho_binned = np.full_like(centers, np.nan)

    for i in range(len(centers)):
        mask = (r >= bins[i]) & (r < bins[i + 1])
        if np.any(mask):
            rho_binned[i] = np.median(rho_sph[mask])

    return centers, rho_binned


def nfw_density(r, M200, c):
    G = 4.30091e-6  # kpc (km/s)^2 / Msun
    H0 = 70.0 / 1e3  # km/s/kpc

    rho_crit = 3 * H0**2 / (8 * np.pi * G)
    R200 = (3 * M200 / (4 * np.pi * 200 * rho_crit)) ** (1 / 3)
    rs = R200 / c

    f = np.log(1 + c) - c / (1 + c)
    rho0 = M200 / (4 * np.pi * rs**3 * f)

    x = r / rs
    return rho0 / (x * (1 + x) ** 2)


def main():
    fig, (ax, ax_res) = plt.subplots(
        2,
        1,
        figsize=(7.2, 7.5),
        sharex=True,
        gridspec_kw={"height_ratios": [3, 1], "hspace": 0},
    )
    lw = 2.0
    ec = "k"
    path_effects = [mpe.Stroke(linewidth=1.8 * lw, foreground=ec), mpe.Normal()]

    with h5py.File(fname, "r") as f:
        pos = np.array(f["SIDMParticles"]["Coordinates"])
        mass = np.array(f["SIDMParticles"]["Masses"]) * 1e10
        rho_sph = np.array(f["SIDMParticles"]["Densities"]) * get_unit_conversions(f)

        time = (
            f["Header"].attrs["Time"][0]
            * f["Units"].attrs["Unit time in cgs (U_t)"][0]
            / 3.15576e16
        )

    r = compute_radius(pos)

    # Binned density
    rs, dens = binned_density(r, mass, rsp)

    # SPH density
    centers, rho_sph_binned = sph_density_profile(r, rho_sph, rsp)

    mask = dens > 0
    rs, dens = rs[mask], dens[mask]

    # Relative residuals w.r.t. the analytic NFW profile the ICs were drawn from
    sph_mask = np.isfinite(rho_sph_binned)
    res_sph = rho_sph_binned[sph_mask] / nfw_density(centers[sph_mask], 5e9, 12) - 1
    res_bin = dens / nfw_density(rs, 5e9, 12) - 1

    # Plot
    ax.plot(
        centers,
        rho_sph_binned,
        lw=lw,
        path_effects=path_effects,
        label=rf"$t={time:.3f}$ Gyr (SPH)",
    )
    ax.plot(rs, dens, lw=lw, path_effects=path_effects, label=rf"$t={time:.3f}$ Gyr")
    ax.plot(
        rsp,
        nfw_density(rsp, 5e9, 12),
        lw=1.2,
        color="black",
        ls="--",
        zorder=0,
        label="NFW",
    )
    ax.set_ylabel(r"$\rho(r)$ [$M_{\odot}$ kpc$^{-3}$]")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.legend()

    # Residual panel
    ax_res.axhline(0.0, color="black", lw=1.2, ls="--", zorder=0)
    ax_res.plot(
        centers[sph_mask], res_sph, lw=lw, color="C0", path_effects=path_effects
    )
    ax_res.plot(rs, res_bin, lw=lw, color="C1", path_effects=path_effects)
    ax_res.set_xlabel("r [kpc]")
    ax_res.set_ylabel(r"$\rho/\rho_{\rm NFW} - 1$")
    ax_res.set_xscale("log")
    ax_res.set_ylim(-0.5, 0.5)
    ax_res.set_yticks([-0.25, 0.0, 0.25])

    plt.tight_layout()
    plt.savefig("density_comparison.png")


if __name__ == "__main__":
    main()
