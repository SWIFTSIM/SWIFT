################################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
Coarse direction/magnitude sanity check for Grackle's photoelectric heating
response to the LW/FUV field this codebase injects: compares a
`with_photoelectric_heating: 1` run against an otherwise-identical
`with_photoelectric_heating: 0` run (same IC, same density/metallicity, same
`time_end`). Bins gas by radius from the star and checks, in every radial
bin where the injected field is non-negligible, that (1) the heating-on run
is hotter than the heating-off run (sign check) and (2) the heating-on
temperature stays under a physically-motivated ceiling rather than showing
signs of a numerical blow-up (magnitude check).

This is deliberately NOT a quantitative match to any PDR literature target
(that is Tier 2's job) -- just "didn't break the sign or blow up".
"""

import argparse
import glob
import sys

import h5py
import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

K_B_CGS = 1.380649e-16
M_P_CGS = 1.67262192e-24
GAMMA_M1 = 2.0 / 3.0  # monatomic ideal gas, gamma = 5/3

# Neither run is expected to sustain a temperature above this: known
# photoelectric-/collisional-ionization-equilibrium temperatures for
# irradiated atomic gas top out around 1e4-1e5 K (Lyman-alpha and
# collisional-ionization cooling both strengthen sharply above ~1e4 K).
# 1e6 K would require a channel this setup does not have and is a much
# more plausible signature of a numerical blow-up than of real heating.
IMPLAUSIBLE_TEMPERATURE_K = 1e6


def parse_options():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--on-snapshot",
        default="snap_on/snapshot_*.hdf5",
        help="Glob pattern for the with_photoelectric_heating=1 run's "
        "snapshots (default: %(default)s)",
    )
    parser.add_argument(
        "--off-snapshot",
        default="snap_off/snapshot_*.hdf5",
        help="Glob pattern for the with_photoelectric_heating=0 run's "
        "snapshots (default: %(default)s)",
    )
    parser.add_argument(
        "--n-bins", type=int, default=12, help="Number of radial bins."
    )
    parser.add_argument(
        "--field-threshold",
        type=float,
        default=1e-3,
        help="A radial bin is considered field-illuminated if the on-run's "
        "mean (u_FUV+u_LW) there exceeds this fraction of the innermost "
        "bin's value (default: %(default)s).",
    )
    parser.add_argument(
        "--output",
        default="isrf_photoelectric_heating_check.png",
        help="Output plot filename.",
    )
    return parser.parse_args()


def load_snapshot(path):
    with h5py.File(path, "r") as f:
        header = f["/Header"]
        time = float(np.asarray(header.attrs["Time"]).flat[0])
        boxsize = np.asarray(header.attrs["BoxSize"], dtype=float).flatten()[0]
        units = f["/Units"]
        unit_length_cgs = float(
            np.asarray(units.attrs["Unit length in cgs (U_L)"]).flat[0]
        )
        unit_time_cgs = float(
            np.asarray(units.attrs["Unit time in cgs (U_t)"]).flat[0]
        )
        u_to_cgs = (unit_length_cgs / unit_time_cgs) ** 2

        gas = f["/PartType0"]
        pos = gas["Coordinates"][:, :]
        u = gas["InternalEnergies"][:] * u_to_cgs
        u_fuv = gas["FUVSpecificEnergies"][:]
        u_lw = gas["LWSpecificEnergies"][:]

        # Species-based mean molecular weight when available
        # (COOLING_GRACKLE_MODE >= 1); otherwise fall back to the
        # primordial-neutral value (X=0.76, Y=0.24, 1/mu = X + Y/4).
        if "HI" in gas:
            inv_mu = (
                gas["HI"][:]
                + 2.0 * gas["HII"][:]
                + 0.25 * gas["HeI"][:]
                + 0.5 * gas["HeII"][:]
                + 0.75 * gas["HeIII"][:]
            )
            mu = 1.0 / np.clip(inv_mu, 1e-10, None)
        else:
            mu = np.full(pos.shape[0], 1.0 / 0.82)

        temperature = u * GAMMA_M1 * mu * M_P_CGS / K_B_CGS

        star = f["/PartType4"]
        star_pos = star["Coordinates"][0, :]

    return dict(
        time=time,
        boxsize=boxsize,
        pos=pos,
        u_fuv=u_fuv,
        u_lw=u_lw,
        temperature=temperature,
        star_pos=star_pos,
    )


def radial_distance(pos, star_pos, boxsize):
    """Minimum-image radial distance from the star, for a periodic box."""
    dx = pos - star_pos
    dx -= boxsize * np.round(dx / boxsize)
    return np.sqrt(np.sum(dx**2, axis=1))


def bin_by_radius(r, values, edges):
    n_bins = len(edges) - 1
    binned = np.full(n_bins, np.nan)
    for i in range(n_bins):
        sel = (r >= edges[i]) & (r < edges[i + 1])
        if sel.sum() > 0:
            binned[i] = np.mean(values[sel])
    return binned


def main():
    opt = parse_options()

    on_files = sorted(glob.glob(opt.on_snapshot))
    off_files = sorted(glob.glob(opt.off_snapshot))
    if not on_files:
        raise RuntimeError(f"No snapshots match {opt.on_snapshot!r}")
    if not off_files:
        raise RuntimeError(f"No snapshots match {opt.off_snapshot!r}")

    # Use the last snapshot of each run: both are designed to reach thermal
    # equilibrium well before time_end (see README).
    on = load_snapshot(on_files[-1])
    off = load_snapshot(off_files[-1])

    r_on = radial_distance(on["pos"], on["star_pos"], on["boxsize"])
    r_off = radial_distance(off["pos"], off["star_pos"], off["boxsize"])

    # Common radial binning, out to the extent of whichever run's
    # illuminated gas has spread furthest (the heating-on run's near-star
    # gas mildly expands as it heats -- see README).
    r_max = max(r_on.max(), r_off.max()) * 0.5
    edges = np.linspace(0.0, r_max, opt.n_bins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])

    field_on = bin_by_radius(r_on, on["u_fuv"] + on["u_lw"], edges)
    T_on = bin_by_radius(r_on, on["temperature"], edges)
    T_off = bin_by_radius(r_off, off["temperature"], edges)

    # A bin counts as "illuminated" if the on-run's own field there is a
    # non-negligible fraction of the peak (innermost, non-nan) bin value --
    # avoids comparing temperatures out where neither run's gas is affected
    # by the star at all (both would trivially sit at the same floor).
    valid = ~np.isnan(field_on) & ~np.isnan(T_on) & ~np.isnan(T_off)
    peak_field = np.nanmax(field_on[valid]) if valid.any() else 0.0
    illuminated = valid & (field_on > opt.field_threshold * peak_field)

    print(f"On-run snapshot: {on_files[-1]} (t={on['time']:.4e})")
    print(f"Off-run snapshot: {off_files[-1]} (t={off['time']:.4e})")
    print(f"{illuminated.sum()}/{opt.n_bins} radial bins are illuminated "
          f"(field > {opt.field_threshold:.1e} x peak).")

    ok_sign = True
    ok_magnitude = True
    for i in np.where(illuminated)[0]:
        sign_ok = T_on[i] > T_off[i]
        magnitude_ok = T_on[i] < IMPLAUSIBLE_TEMPERATURE_K
        ok_sign &= sign_ok
        ok_magnitude &= magnitude_ok
        print(
            f"  bin {i}: r={centres[i]:.4e}, T_on={T_on[i]:.4e} K, "
            f"T_off={T_off[i]:.4e} K, "
            f"{'PASS' if sign_ok else 'FAIL'} (sign), "
            f"{'PASS' if magnitude_ok else 'FAIL'} (magnitude)"
        )

    if illuminated.sum() == 0:
        print("No illuminated bins found -- cannot check anything.")
        ok_sign = ok_magnitude = False

    print(f"Sign check (T_on > T_off in every illuminated bin): "
          f"{'PASS' if ok_sign else 'FAIL'}")
    print(f"Magnitude check (T_on < {IMPLAUSIBLE_TEMPERATURE_K:.1e} K in "
          f"every illuminated bin): {'PASS' if ok_magnitude else 'FAIL'}")

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.semilogy(centres[valid], T_on[valid], "o-", label="heating on")
    ax.semilogy(centres[valid], T_off[valid], "s-", label="heating off")
    ax.axvspan(
        centres[illuminated].min() if illuminated.any() else 0,
        centres[illuminated].max() if illuminated.any() else 0,
        color="grey",
        alpha=0.15,
        label="illuminated bins",
    )
    ax.set_xlabel("r (internal length units)")
    ax.set_ylabel("Temperature (K)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(opt.output, dpi=150)
    print(f"Plot saved to {opt.output}")

    if not (ok_sign and ok_magnitude):
        sys.exit(1)


if __name__ == "__main__":
    main()
