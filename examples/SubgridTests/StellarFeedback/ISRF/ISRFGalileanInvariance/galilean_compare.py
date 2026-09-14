#!/usr/bin/env python3
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
Galilean-invariance check of the ISRF propagation scheme.

The scheme evolves the mass-specific FUV/LW energy u and flux F along the
gas trajectories and reads no particle velocity, so F is the fluid-frame
flux. A run whose gas and star share a uniform bulk velocity V must then
carry, on every gas particle, the same u and F as the run at rest, up to
round-off: the boosted particle i is the same material element as the rest
particle i, displaced by V t. Particles are matched by ParticleIDs.

Metrics, per band, maximum over every common snapshot:
  D_u = max_i |u_boost - u_rest| / max_i |u_rest|
  D_F = max_i |F_boost - F_rest| / max_i |F_rest|, over particles with
        |F_rest| > F_REL_FLOOR * max|F_rest| (the excluded count is printed).
The noise floor N0(band, metric) is the larger of the same metric for two
controls at rest: an identical repeat, and the same run with every position
shifted by a constant (different cell layout, same physics).

PASS: every boosted run has D <= NOISE_FACTOR * N0 for each band and metric.
FAIL: otherwise, or any non-finite value in a compared field.
D_u is also printed against U_REFERENCE, the round-off level of a rigid boost
without periodic wrap.
"""

import argparse
import glob
import json
import os
import sys

import h5py
import numpy as np

NOISE_FACTOR = 3.0
U_REFERENCE = 1e-5
F_REL_FLOOR = 1e-6
BANDS = ("FUV", "LW")


def parse_options():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--measure-c-hyp",
        type=str,
        default=None,
        help="Print the median closure c_hyp (km/s) of this run and exit.",
    )
    parser.add_argument("--rest", type=str, help="Run directory at rest.")
    parser.add_argument("--repeat", type=str, help="Identical repeat of the rest run.")
    parser.add_argument(
        "--shifted", type=str, help="Rest run with a uniform position shift."
    )
    parser.add_argument(
        "--boosted",
        type=str,
        nargs="*",
        default=[],
        help="Run directories with a uniform bulk velocity.",
    )
    parser.add_argument("--json-out", type=str, default="galilean_summary.json")
    return parser.parse_args()


def load_used_parameters(run):
    import yaml

    with open(os.path.join(run, "used_parameters.yml")) as f:
        return yaml.safe_load(f)


def modal_bulk_dt(run):
    """Most frequent interval between successive steps in timesteps.txt."""
    times = []
    with open(os.path.join(run, "timesteps.txt")) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            if int(parts[0]) == 0:
                continue
            times.append(float(parts[1]))
    diffs = np.round(np.diff(np.array(times)), 14)
    vals, counts = np.unique(diffs, return_counts=True)
    return float(vals[np.argmax(counts)])


def snapshots(run):
    return sorted(glob.glob(os.path.join(run, "snap", "snapshot_*.hdf5")))


def load_snapshot(path):
    with h5py.File(path, "r") as f:
        gas = f["/PartType0"]
        order = np.argsort(gas["ParticleIDs"][:])
        snap = dict(
            time=float(np.asarray(f["/Header"].attrs["Time"]).flat[0]),
            boxsize=float(np.asarray(f["/Header"].attrs["BoxSize"]).flat[0]),
            ids=gas["ParticleIDs"][:][order],
            pos=gas["Coordinates"][:][order].astype(np.float64),
            h=gas["SmoothingLengths"][:][order].astype(np.float64),
        )
        for band in BANDS:
            snap[f"u_{band}"] = gas[f"{band}SpecificEnergies"][:][order].astype(
                np.float64
            )
            snap[f"F_{band}"] = gas[f"{band}SpecificFluxes"][:][order].astype(
                np.float64
            )
    return snap


def measure_c_hyp(run):
    params = load_used_parameters(run)
    margin = float(params["GEARFeedback"].get("ISRF_c_hyp_margin", 0.5))
    h_med = float(np.median(load_snapshot(snapshots(run)[-1])["h"]))
    return margin * h_med / modal_bulk_dt(run)


def bulk_velocity(run):
    with h5py.File(os.path.join(run, "ICs_isrf_galilean_invariance.hdf5"), "r") as f:
        return np.asarray(f["/PartType0/Velocities"][:]).astype(np.float64).mean(axis=0)


def compare(run_a, run_b, velocity):
    """Return per-band D_u, D_F, the position residual in h and a finiteness flag."""
    files_a, files_b = snapshots(run_a), snapshots(run_b)
    if len(files_a) != len(files_b) or len(files_a) == 0:
        raise RuntimeError(
            f"{run_b}: {len(files_b)} snapshots, rest has {len(files_a)}"
        )
    out = {
        band: dict(D_u=0.0, D_F=0.0, n_F_excluded=0, min_abs_F=np.inf) for band in BANDS
    }
    finite = True
    max_pos_residual_h = 0.0
    for fa, fb in zip(files_a, files_b):
        a, b = load_snapshot(fa), load_snapshot(fb)
        if not np.array_equal(a["ids"], b["ids"]):
            raise RuntimeError(f"{fb}: particle IDs differ from {fa}")
        if abs(a["time"] - b["time"]) > 1e-9 * max(abs(a["time"]), 1e-30):
            raise RuntimeError(f"{fb}: time {b['time']} differs from {a['time']}")
        L = a["boxsize"]
        dx = b["pos"] - a["pos"] - velocity[None, :] * a["time"]
        dx -= L * np.round(dx / L)
        max_pos_residual_h = max(
            max_pos_residual_h,
            float(np.max(np.sqrt(np.sum(dx**2, axis=1)) / np.median(a["h"]))),
        )
        for band in BANDS:
            ua, ub = a[f"u_{band}"], b[f"u_{band}"]
            Fa, Fb = a[f"F_{band}"], b[f"F_{band}"]
            for arr in (ua, ub, Fa, Fb):
                finite &= bool(np.all(np.isfinite(arr)))
            if not finite:
                continue
            u_scale = np.max(np.abs(ua))
            if u_scale > 0.0:
                out[band]["D_u"] = max(
                    out[band]["D_u"], float(np.max(np.abs(ub - ua)) / u_scale)
                )
            Fa_norm = np.sqrt(np.sum(Fa**2, axis=1))
            F_scale = np.max(Fa_norm)
            if F_scale > 0.0:
                keep = Fa_norm > F_REL_FLOOR * F_scale
                dF = np.sqrt(np.sum((Fb - Fa) ** 2, axis=1))
                out[band]["D_F"] = max(
                    out[band]["D_F"], float(np.max(dF[keep]) / F_scale)
                )
                out[band]["n_F_excluded"] = max(
                    out[band]["n_F_excluded"], int(np.sum(~keep))
                )
                nonzero = Fa_norm[Fa_norm > 0.0]
                if nonzero.size:
                    out[band]["min_abs_F"] = min(
                        out[band]["min_abs_F"], float(nonzero.min())
                    )
    return out, max_pos_residual_h, finite


def main():
    opt = parse_options()
    if opt.measure_c_hyp is not None:
        print(f"{measure_c_hyp(opt.measure_c_hyp):.6g}")
        return

    c_hyp = measure_c_hyp(opt.rest)
    print(f"Median closure c_hyp (rest run): {c_hyp:.4f} km/s")
    summary = dict(c_hyp_kms=c_hyp, noise_factor=NOISE_FACTOR, runs={}, noise_floor={})
    all_finite = True

    N0 = {band: dict(D_u=0.0, D_F=0.0) for band in BANDS}
    for control in (opt.repeat, opt.shifted):
        noise, _, finite = compare(opt.rest, control, np.zeros(3))
        all_finite &= finite
        for band in BANDS:
            for metric in ("D_u", "D_F"):
                N0[band][metric] = max(N0[band][metric], noise[band][metric])
            print(
                f"noise floor ({control}) {band}: D_u={noise[band]['D_u']:.3e} "
                f"D_F={noise[band]['D_F']:.3e}"
            )
        summary["noise_floor"][control] = dict(bands=noise, finite=finite)
    summary["N0"] = N0
    for band in BANDS:
        print(
            f"bar {band}: D_u <= {NOISE_FACTOR * N0[band]['D_u']:.3e}, "
            f"D_F <= {NOISE_FACTOR * N0[band]['D_F']:.3e}"
        )

    passed = all_finite
    for run in opt.boosted:
        velocity = bulk_velocity(run)
        res, pos_res_h, finite = compare(opt.rest, run, velocity)
        all_finite &= finite
        run_pass = finite
        for band in BANDS:
            ok = all(
                res[band][metric] <= NOISE_FACTOR * N0[band][metric]
                for metric in ("D_u", "D_F")
            )
            run_pass &= ok
            print(
                f"{run} V={velocity[0]:.4g} km/s ({velocity[0] / c_hyp:.3g} c_hyp) {band}: "
                f"D_u={res[band]['D_u']:.3e} "
                f"({'<=' if res[band]['D_u'] <= U_REFERENCE else '>'} {U_REFERENCE:.0e}) "
                f"D_F={res[band]['D_F']:.3e} "
                f"(F excluded {res[band]['n_F_excluded']}, min|F| {res[band]['min_abs_F']:.3e}) "
                f"{'ok' if ok else 'ABOVE BAR'}"
            )
        print(
            f"{run}: max position residual |x_b - x_a - V t| = {pos_res_h:.3e} h, "
            f"finite={finite}"
        )
        passed &= run_pass
        summary["runs"][run] = dict(
            velocity_kms=velocity.tolist(),
            bands=res,
            position_residual_h=pos_res_h,
            finite=finite,
            passed=bool(run_pass),
        )

    if not all_finite:
        print("Non-finite value in a compared field.")
    verdict = "PASS" if passed and all_finite else "FAIL"
    summary["verdict"] = verdict
    with open(opt.json_out, "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print(f"Galilean invariance: {verdict}")
    sys.exit(0 if verdict == "PASS" else 1)


if __name__ == "__main__":
    main()
