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
Tier 1 ("identity") check: compares the cosmological run's simulated
r_hii(t) directly against the non-cosmological ../StromgrenSphere control's
own r_hii(t) -- not against the analytic formula (see
cosmo_stromgren_analytic_check.py for that, used for Tier 2). A cosmological
frame bug (a missing/spurious scale-factor factor, or a d(ln a)-vs-proper-time
confusion) shows up here as a deviation from the control even though the
underlying physics is unchanged; the analytic formula is not needed to catch
that.

Both runs are read through read_simulated_r_hii() from
cosmo_stromgren_analytic_check.py, which always converts comoving snapshot
fields to physical using each field's own a-scale exponent -- for the
non-cosmological control this is a no-op (a=1 throughout), so the same code
path handles both runs identically.

Usage:
    python3 cosmo_identity_check.py --control control/snap/snapshot_*.hdf5 \\
        --cosmo identity/snap/snapshot_*.hdf5
"""

import argparse
import glob

import h5py

import numpy as np
from astropy import units as u

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cosmo_stromgren_analytic_check import read_simulated_r_hii


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--control",
        required=True,
        help="Glob pattern for the non-cosmological control run's snapshots.",
    )
    parser.add_argument(
        "--cosmo",
        required=True,
        help="Glob pattern for the cosmological (tier=identity) run's snapshots.",
    )
    parser.add_argument(
        "-o",
        "--output",
        default="cosmo_identity_check.png",
        help="Output plot filename.",
    )
    parser.add_argument(
        "--tol",
        type=float,
        default=0.1,
        help="Relative-error tolerance for a pass (default: 0.1). Generous "
        "enough to absorb the two runs' independently-sampled particle "
        "counts/positions (same RNG seed, but a slightly different box "
        "size and particle count -- see makeIC.py), while still catching "
        "an O(1) frame bug.",
    )
    args = parser.parse_args()

    if not glob.glob(args.control):
        raise RuntimeError(f"No control snapshots found matching {args.control!r}.")
    if not glob.glob(args.cosmo):
        raise RuntimeError(f"No cosmological snapshots found matching {args.cosmo!r}.")

    t_ctrl, r_ctrl, m_ctrl, n_H_ctrl = read_simulated_r_hii(args.control)
    t_cosmo, r_cosmo, m_cosmo, n_H_cosmo = read_simulated_r_hii(args.cosmo)

    print(
        f"Control : star={m_ctrl:.3f} Msun, n_H(t=0)={n_H_ctrl:.4g}, "
        f"{len(t_ctrl)} snapshots, t in [{t_ctrl.min():.4g}, {t_ctrl.max():.4g}]"
    )
    print(
        f"Cosmo   : star={m_cosmo:.3f} Msun, n_H(t=0)={n_H_cosmo:.4g}, "
        f"{len(t_cosmo)} snapshots, t in [{t_cosmo.min():.4g}, {t_cosmo.max():.4g}]"
    )

    # Compare at the (coarser) cosmological run's own sample times, linearly
    # interpolating the control run's finer-sampled r_hii(t) onto them.
    # Only the overlapping, both-alive (r>0) time range is meaningful.
    t_lo = max(t_ctrl.min(), t_cosmo.min())
    t_hi = min(t_ctrl.max(), t_cosmo.max())
    in_range = (t_cosmo >= t_lo) & (t_cosmo <= t_hi) & (r_cosmo > 0)
    if not np.any(in_range):
        raise RuntimeError(
            "No overlapping, alive-star time range between the control and "
            "cosmological runs -- check both reached a comparable number of "
            "HII rebuild passes before comparing."
        )
    r_ctrl_interp = (
        np.interp(
            t_cosmo.to(u.Myr).value[in_range],
            t_ctrl.to(u.Myr).value,
            r_ctrl.to(u.kpc).value,
        )
        * u.kpc
    )
    r_cosmo_valid = r_cosmo[in_range]

    # Two regimes carry no frame information and must not enter the verdict.
    # Saturation: this setup's equilibrium Stromgren radius exceeds the box
    # half-width, so at late times BOTH runs tag the whole box and the
    # max-tagged-radius pins at the periodic corner distance -- comparing
    # there measures corner-pocket tag churn, not cosmology. Shot noise: at
    # radii of a few interparticle spacings, the max-radius estimator jumps
    # by whole spacings, so relative errors there are discretization noise.
    with h5py.File(sorted(glob.glob(args.control))[0], "r") as h:
        box_kpc = float(h["Header"].attrs["BoxSize"][0])
        n_gas_1d = float(h["Header"].attrs["NumPart_Total"][0]) ** (1.0 / 3.0)
    r_saturation = 0.35 * box_kpc * u.kpc
    r_floor = 2.0 * (box_kpc / n_gas_1d) * u.kpc
    also_alive = (
        (r_ctrl_interp > r_floor)
        & (r_ctrl_interp < r_saturation)
        & (r_cosmo_valid < r_saturation)
    )
    n_sat = int(np.sum(r_ctrl_interp >= r_saturation))
    n_tiny = int(np.sum((r_ctrl_interp > 0) & (r_ctrl_interp <= r_floor)))
    if not np.any(also_alive):
        raise RuntimeError(
            "No comparison points survive the saturation/shot-noise cuts -- "
            "the growth phase was not sampled; increase the snapshot cadence."
        )
    rel_error = (
        np.abs(r_cosmo_valid[also_alive] - r_ctrl_interp[also_alive])
        / r_ctrl_interp[also_alive]
    )
    rel_error = rel_error.decompose().value

    max_err = float(np.max(rel_error))
    mean_err = float(np.mean(rel_error))
    verdict = "PASS" if mean_err <= args.tol else "FAIL"
    print(
        f"r_hii(t) relative error over {int(np.sum(also_alive))} growth-phase "
        f"points (excluded: {n_sat} saturated, {n_tiny} below the "
        f"shot-noise floor): mean={mean_err:.2%}  max={max_err:.2%}  "
        f"(tol={args.tol:.0%} on the mean)  [{verdict}]"
    )

    fig, ax = plt.subplots()
    ax.plot(t_ctrl.to(u.Myr), r_ctrl.to(u.pc), "-", label="Control (non-cosmological)")
    ax.plot(
        t_cosmo.to(u.Myr), r_cosmo.to(u.pc), "o", label="Cosmological (tier=identity)"
    )
    ax.set_xlabel("Time since a_begin [Myr]")
    ax.set_ylabel("HII region radius [pc] (physical)")
    ax.legend()
    ax.grid(True, linestyle="--", alpha=0.5)
    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"Plot saved to {args.output}")


if __name__ == "__main__":
    main()
