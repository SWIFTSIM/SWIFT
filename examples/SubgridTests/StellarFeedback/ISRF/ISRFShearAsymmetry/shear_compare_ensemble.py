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
Ensemble, time-averaged gate for ISRFShearAsymmetry: gates a merge on the
A_centroid/A_spread/A_energy reflection-symmetry metrics (`--mode gate` in
shear_compare.py remains for a quick single-run/single-control check).
Reuses isrf_shear_asymmetry_check.py's per-band moment functions, so the
metric definition itself has one source. Blob geometry, variant=shear only
(the two moving blobs whose reflection asymmetry the gate probes); slab
geometry is out of scope.

Design:
(a) each realisation's asymmetry metric is TIME-AVERAGED over the settled
    half of the run, t >= WINDOW_FRAC * t_last (default 0.5): burn in half
    the run, generic and not tuned to any one fixture's outcome.
(b) the bar comes from an ENSEMBLE of N zero-shear control realisations
    (default 4), each on an independently glass-shifted IC (makeIC.py
    --glass-shift), not from one control run's single value.
(c) the signal is the mean over N sheared (v_shear > 0) realisations (same
    glass-shift protocol). A v_shear < 0 realisation is NOT pooled into
    this group and is NOT sign-flipped: A_centroid is a SUM
    (Dx_A + Dx_B)/h, not a difference, so it does not change sign under an
    A<->B blob relabelling the way A_spread/A_energy do, and a v -> -v run
    is not generally the same physical configuration reflected (it depends
    on which symmetry maps +v to -v: a label swap by periodic (L/2, L/2)
    translation and an x-mirror have different parities on these metrics).
    `--minus-runs` accepts v_shear < 0 realisations and reports them as
    their own group (mean, SE) with a separate, report-only comparison
    against the control ensemble; they never enter the shear/control gate.
(d) gate: |mean(shear) - mean(control)| <= 3*SE_diff,
    SE_diff = sqrt(SE(shear)^2 + SE(control)^2), SE = sample std (ddof=1)
    / sqrt(N). With N < 2 realisations in a group the standard error
    cannot be estimated: that band/metric is reported INCONCLUSIVE, never
    silently PASS or FAIL. A non-finite time-average in any realisation
    FAILs that band/metric outright (the offending run is named in the
    output), never silently dropped from the group.

VOID handling (KH contamination at the last snapshot, negative-weight
share) is unchanged from isrf_shear_asymmetry_check.py: a void realisation
is dropped from its group and reported, not averaged in.
"""

import argparse
import glob
import json
import os
import sys

import h5py
import numpy as np

from isrf_shear_asymmetry_check import (
    asymmetry_from_pair,
    blob_membership,
    blob_pair_moments,
    kh_contamination_void,
)

WINDOW_FRAC = 0.5
SE_FACTOR = 3.0
# v_shear is measured from particle velocities (see realisation_time_series),
# so a "zero-shear" control does not read exactly 0; below this threshold
# it is treated as zero shear (well under the default 1 km/s v_shear and
# far above the ~1e-5 km/s measurement noise seen on this fixture).
SHEAR_DETECT_KMS = 1e-3


def parse_options():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--shear-runs",
        nargs="+",
        required=True,
        help="Sheared (v_shear > 0) realisation run dirs.",
    )
    parser.add_argument(
        "--minus-runs",
        nargs="+",
        default=[],
        help="v_shear < 0 realisation run dirs, reported as their own group "
        "(never pooled into --shear-runs; see the module docstring).",
    )
    parser.add_argument(
        "--control-runs",
        nargs="+",
        required=True,
        help="Zero-shear control realisation run dirs (independent glass shifts).",
    )
    parser.add_argument("--window-frac", type=float, default=WINDOW_FRAC)
    parser.add_argument(
        "--pulse-sigma-h",
        type=float,
        default=2.0,
        help="Must match the value makeIC.py/isrf_shear_asymmetry_check.py used.",
    )
    parser.add_argument(
        "--layer-width-h",
        type=float,
        default=4.0,
        help="Must match the value makeIC.py's --layer-width-h used (its own default).",
    )
    parser.add_argument("--json-out", default="shear_ensemble_metrics.json")
    return parser.parse_args()


def load_snapshot(path):
    with h5py.File(path, "r") as f:
        time = float(np.asarray(f["/Header"].attrs["Time"]).flat[0])
        boxsize = np.asarray(f["/Header"].attrs["BoxSize"], dtype=float).flatten()[0]
        gas = f["/PartType0"]
        pos = gas["Coordinates"][:, :]
        vel = gas["Velocities"][:, :]
        mass = gas["Masses"][:].astype(np.float64)
        rho = gas["Densities"][:].astype(np.float64)
        h = gas["SmoothingLengths"][:].astype(np.float64)
        ids = gas["ParticleIDs"][:]
        u_fuv = gas["FUVSpecificEnergies"][:].astype(np.float64)
        u_lw = gas["LWSpecificEnergies"][:].astype(np.float64)
    return dict(
        time=time,
        boxsize=boxsize,
        pos=pos,
        vel=vel,
        mass=mass,
        rho=rho,
        h=h,
        ids=ids,
        u_fuv=u_fuv,
        u_lw=u_lw,
    )


def realisation_time_series(run_dir, pulse_sigma_h, layer_width_h):
    """Per-snapshot A_centroid/A_spread/A_energy (both bands), the run's own
    v_shear, void flag and h_med_last, for one realisation directory.
    `pulse_sigma_h`/`layer_width_h` must match the values used to generate
    the run's IC (makeIC.py's own --pulse-sigma-h/--layer-width-h)."""
    files = sorted(glob.glob(os.path.join(run_dir, "snap", "snapshot_*.hdf5")))
    if not files:
        raise RuntimeError(f"No snapshots in {run_dir}")

    snap0 = load_snapshot(files[0])
    snap_last = load_snapshot(files[-1])
    order0 = np.argsort(snap0["ids"])
    pos0 = snap0["pos"][order0]
    mass0 = snap0["mass"][order0]
    vel0 = snap0["vel"][order0]
    L = snap0["boxsize"]
    n_gas = pos0.shape[0]
    h_mean_analytic = 1.2348 * L / n_gas ** (1.0 / 3.0)
    sigma = pulse_sigma_h * h_mean_analytic
    in_A, in_B, centre_A, centre_B = blob_membership(pos0, L, sigma)
    # v_shear is not in used_parameters.yml (it is a makeIC.py IC-generation
    # option, not a runtime SWIFT parameter): recover it by a linear fit of
    # the t=0 vx(y) to makeIC.py's own tanh profile (fits the whole profile
    # rather than differencing blob averages, which the unsaturated tanh
    # tails would bias low).
    y_frac = pos0[:, 1] / L
    d_frac = layer_width_h * h_mean_analytic / L
    basis = 0.5 * (
        np.tanh((y_frac - 0.25) / d_frac) - np.tanh((y_frac - 0.75) / d_frac) - 1.0
    )
    v_shear = float(np.sum(vel0[:, 0] * basis) / np.sum(basis**2))

    orderL = np.argsort(snap_last["ids"])
    if not np.array_equal(snap0["ids"][order0], snap_last["ids"][orderL]):
        raise RuntimeError(f"{run_dir}: ParticleIDs mismatch first vs last snapshot.")
    h_med_last = float(np.median(snap_last["h"]))

    # KH-contamination / density-drift validity, last snapshot only (this
    # check is inherently instantaneous; unrelated to the time-averaging
    # redesign, kept as isrf_shear_asymmetry_check.py defines it).
    rho0 = snap0["rho"][order0]
    rho_last = snap_last["rho"][orderL]
    is_sheared = abs(v_shear) > SHEAR_DETECT_KMS
    vy_rms = float(np.sqrt(np.mean(snap_last["vel"][:, 1] ** 2)))
    kh = vy_rms / abs(v_shear) if is_sheared else float("nan")
    drho = float(np.max(np.abs(rho_last - rho0) / rho0))
    kh_void = kh_contamination_void(kh, drho, is_sheared)

    times, series = [], {
        b: {"A_centroid": [], "A_spread": [], "A_energy": []} for b in ("FUV", "LW")
    }
    neg_weight_void = False
    for fn in files:
        s = load_snapshot(fn)
        order = np.argsort(s["ids"])
        pos = s["pos"][order]
        t = s["time"]
        times.append(t)
        for band, key in (("FUV", "u_fuv"), ("LW", "u_lw")):
            u = s[key][order]
            pair = blob_pair_moments(
                pos, u, mass0, in_A, in_B, centre_A, centre_B, v_shear, t, L
            )
            asym = asymmetry_from_pair(pair, h_med_last)
            neg_weight_void |= asym["void_neg_weight"]
            series[band]["A_centroid"].append(asym["A_centroid"])
            series[band]["A_spread"].append(asym["A_spread"])
            series[band]["A_energy"].append(asym["A_energy"])

    void = bool(kh_void or neg_weight_void)
    return dict(
        run_dir=run_dir,
        times=np.array(times),
        series=series,
        v_shear=v_shear,
        void=void,
        void_reason=(
            "KH/drho" if kh_void else ("neg-weight" if neg_weight_void else None)
        ),
    )


def time_average(times, values, window_frac):
    times = np.asarray(times)
    values = np.asarray(values, dtype=float)
    t_last = times[-1]
    sel = times >= window_frac * t_last
    if not np.any(sel):
        sel = np.zeros(len(times), dtype=bool)
        sel[-1] = True
    vals = values[sel]
    n_used = int(sel.sum())
    if not np.all(np.isfinite(vals)):
        return np.nan, n_used
    return float(np.mean(vals)), n_used


def group_stats(labeled_means):
    """labeled_means: list of (run_dir, time-averaged value). A non-finite
    value is reported in `non_finite_dirs`, never silently averaged out."""
    finite = [(d, v) for d, v in labeled_means if np.isfinite(v)]
    non_finite_dirs = [d for d, v in labeled_means if not np.isfinite(v)]
    vals = np.array([v for _, v in finite])
    n = len(vals)
    mean = float(np.mean(vals)) if n >= 1 else np.nan
    se = float(np.std(vals, ddof=1) / np.sqrt(n)) if n >= 2 else np.nan
    return mean, se, n, non_finite_dirs


def main():
    opt = parse_options()

    shear_series = [
        realisation_time_series(r, opt.pulse_sigma_h, opt.layer_width_h)
        for r in opt.shear_runs
    ]
    minus_series = [
        realisation_time_series(r, opt.pulse_sigma_h, opt.layer_width_h)
        for r in opt.minus_runs
    ]
    control_series = [
        realisation_time_series(r, opt.pulse_sigma_h, opt.layer_width_h)
        for r in opt.control_runs
    ]

    for label, group in (
        ("shear", shear_series),
        ("minus", minus_series),
        ("control", control_series),
    ):
        for r in group:
            status = "VOID" if r["void"] else "ok"
            print(
                f"{label}: {r['run_dir']} v_shear={r['v_shear']:.3g} km/s "
                f"n_snap={len(r['times'])} -> {status}"
                + (f" ({r['void_reason']})" if r["void"] else "")
            )

    results = dict(window_frac=opt.window_frac, bands={})
    overall = "PASS"
    for band in ("FUV", "LW"):
        results["bands"][band] = {}
        print(f"\n=== {band} ===")
        for metric in ("A_centroid", "A_spread", "A_energy"):
            shear_means, minus_means, control_means = [], [], []
            for r in shear_series:
                if r["void"]:
                    continue
                avg, _ = time_average(
                    r["times"], r["series"][band][metric], opt.window_frac
                )
                shear_means.append((r["run_dir"], avg))
            for r in minus_series:
                if r["void"]:
                    continue
                avg, _ = time_average(
                    r["times"], r["series"][band][metric], opt.window_frac
                )
                minus_means.append((r["run_dir"], avg))
            for r in control_series:
                if r["void"]:
                    continue
                avg, _ = time_average(
                    r["times"], r["series"][band][metric], opt.window_frac
                )
                control_means.append((r["run_dir"], avg))

            mean_s, se_s, n_s, bad_s = group_stats(shear_means)
            mean_c, se_c, n_c, bad_c = group_stats(control_means)
            mean_m, se_m, n_m, bad_m = group_stats(minus_means)
            bad_gate = bad_s + bad_c
            if bad_gate:
                # Non-finite realisation feeding the gate: FAIL outright,
                # never average over only the finite ones.
                verdict = "FAIL"
                overall = "FAIL"
                diff = se_diff = float("nan")
                print(f"{metric}: non-finite time-average in {bad_gate} -> FAIL")
            else:
                diff = (
                    abs(mean_s - mean_c)
                    if np.isfinite(mean_s) and np.isfinite(mean_c)
                    else np.nan
                )
                se_diff = (
                    float(np.hypot(se_s, se_c))
                    if np.isfinite(se_s) and np.isfinite(se_c)
                    else np.nan
                )
                if not np.isfinite(diff):
                    verdict = "FAIL"  # non-finite input: never a silent pass
                    overall = "FAIL"
                elif not np.isfinite(se_diff):
                    verdict = "INCONCLUSIVE"
                elif diff <= SE_FACTOR * se_diff:
                    verdict = "PASS"
                else:
                    verdict = "FAIL"
                    overall = "FAIL"
                print(
                    f"{metric}: shear mean={mean_s:+.4e} se={se_s:.4e} (N={n_s})  "
                    f"control mean={mean_c:+.4e} se={se_c:.4e} (N={n_c})  "
                    f"|diff|={diff:.4e}  bar={SE_FACTOR}*se_diff={SE_FACTOR * se_diff if np.isfinite(se_diff) else float('nan'):.4e}"
                    f"  -> {verdict}"
                )
            metric_result = dict(
                mean_shear=mean_s,
                se_shear=se_s,
                n_shear=n_s,
                mean_control=mean_c,
                se_control=se_c,
                n_control=n_c,
                diff=diff,
                se_diff=se_diff,
                verdict=verdict,
            )
            if bad_m:
                print(
                    f"{metric}: [report-only] non-finite minus time-average in {bad_m}"
                )
                metric_result["minus"] = dict(
                    mean=mean_m, se=se_m, n=n_m, non_finite=bad_m
                )
            elif minus_series:
                diff_m = (
                    abs(mean_m - mean_c)
                    if np.isfinite(mean_m) and np.isfinite(mean_c)
                    else np.nan
                )
                se_diff_m = (
                    float(np.hypot(se_m, se_c))
                    if np.isfinite(se_m) and np.isfinite(se_c)
                    else np.nan
                )
                print(
                    f"{metric}: [report-only, not gated] minus mean={mean_m:+.4e} "
                    f"se={se_m:.4e} (N={n_m})  vs control  |diff|={diff_m:.4e}  "
                    f"3*se_diff={SE_FACTOR * se_diff_m if np.isfinite(se_diff_m) else float('nan'):.4e}"
                )
                metric_result["minus"] = dict(
                    mean=mean_m,
                    se=se_m,
                    n=n_m,
                    diff_vs_control=diff_m,
                    se_diff_vs_control=se_diff_m,
                )
            results["bands"][band][metric] = metric_result
    if overall == "PASS" and any(
        results["bands"][b][m]["verdict"] == "INCONCLUSIVE"
        for b in results["bands"]
        for m in results["bands"][b]
    ):
        overall = "INCONCLUSIVE"
    results["overall"] = overall
    with open(opt.json_out, "w") as f:
        json.dump(results, f, indent=2, default=float)
    print(f"\nOverall: {overall}")
    sys.exit(0 if overall == "PASS" else 1)


if __name__ == "__main__":
    main()
