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
Ensemble, time-averaged gate for ISRFShearAsymmetry (redesign of
shear_compare.py's mode=gate; NEW SCRIPT, not a drop-in replacement -- the
two must be reconciled when this branch is merged, see the item 5 follow-up
log). Blob geometry, variant=shear only (the two moving blobs whose
reflection asymmetry the gate probes); slab geometry is out of scope.

Why a redesign (ISRF_MASTER.md 4d, 2026-09-15): the original gate compared
one sheared run's LAST-snapshot A_centroid/A_spread/A_energy against 3x a
single zero-shear control's OWN last-snapshot value. The control's
A_centroid swings between about -0.1 and +0.1 over the run (pure numerical
noise: there is no shear to break the reflection symmetry), so a
single-snapshot bar depends on which snapshot it lands on and which glass
realisation set it: the failing band switched between realisations, and
+v/-v were not mirrors of each other (A_pair 0.55/0.61).

Redesign:
(a) each realisation's asymmetry metric is TIME-AVERAGED over the settled
    half of the run, t >= WINDOW_FRAC * t_last (default 0.5), a rule fixed
    before looking at any A_centroid number: burn in half the run, generic
    and not tuned to this fixture's outcome.
(b) the bar comes from an ENSEMBLE of N zero-shear control realisations
    (default 4), each on an independently glass-shifted IC (makeIC.py
    --glass-shift), not from one control run's single value.
(c) the signal is the mean over N sheared realisations (same glass-shift
    protocol). A realisation with v_shear < 0 is sign-flipped before
    pooling: reflection symmetry v -> -v flips the sign of every odd-in-v
    metric here (A_centroid, A_spread, A_energy all change sign under
    A<->B blob relabelling, which is what v -> -v does), so a -v run
    measures the same physical asymmetry as a +v run on a mirrored glass.
(d) gate: |mean(shear) - mean(control)| <= 3*SE_diff,
    SE_diff = sqrt(SE(shear)^2 + SE(control)^2), SE = sample std (ddof=1)
    / sqrt(N). With N < 2 realisations in a group the standard error
    cannot be estimated: that band/metric is reported INCONCLUSIVE, never
    silently PASS or FAIL.

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
import yaml

NEG_WEIGHT_VOID_THRESHOLD = 0.10
WINDOW_FRAC = 0.5
SE_FACTOR = 3.0
# v_shear is measured from particle velocities (see realisation_time_series),
# so a "zero-shear" control reads a few km/s numerical noise, not exactly 0;
# below this it is treated as zero (well under the default 1 km/s v_shear
# and far above the ~1e-5 km/s measurement noise seen on this fixture).
SHEAR_DETECT_KMS = 1e-3


def parse_options():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--shear-runs", nargs="+", required=True, help="Sheared realisation run dirs."
    )
    parser.add_argument(
        "--control-runs",
        nargs="+",
        required=True,
        help="Zero-shear control realisation run dirs (independent glass shifts).",
    )
    parser.add_argument("--window-frac", type=float, default=WINDOW_FRAC)
    parser.add_argument("--json-out", default="shear_ensemble_metrics.json")
    return parser.parse_args()


def min_image(dx, boxsize):
    return dx - boxsize * np.round(dx / boxsize)


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


def blob_membership(pos0, boxsize, sigma):
    centre_A = np.array([0.25 * boxsize, 0.50 * boxsize, 0.50 * boxsize])
    centre_B = np.array([0.75 * boxsize, 0.00 * boxsize, 0.50 * boxsize])
    dA = min_image(pos0 - centre_A, boxsize)
    dB = min_image(pos0 - centre_B, boxsize)
    in_A = np.sqrt(np.sum(dA**2, axis=1)) <= 3.0 * sigma
    in_B = np.sqrt(np.sum(dB**2, axis=1)) <= 3.0 * sigma
    return in_A, in_B, centre_A, centre_B


def moments(w, dx_x):
    w_sum = w.sum()
    if w_sum == 0:
        return 0.0, 0.0
    return float(np.sum(w * dx_x) / w_sum), float(np.sum(w * dx_x**2) / w_sum)


def realisation_time_series(run_dir):
    """Per-snapshot A_centroid/A_spread/A_energy (both bands), the run's own
    v_shear, void flag and h_med_last, for one realisation directory."""
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
    sigma = 2.0 * h_mean_analytic
    in_A, in_B, centre_A, centre_B = blob_membership(pos0, L, sigma)
    # v_shear is not in used_parameters.yml (it is a makeIC.py IC-generation
    # option, not a runtime SWIFT parameter): recover it by a linear fit of
    # the t=0 vx(y) to makeIC.py's own tanh profile (assumes the default
    # --layer-width-h=4.0; the profile is wide enough at this resolution
    # that a blob-average is biased low by the unsaturated tanh tails, a
    # crude two-bin difference measured ~0.73-0.84 against a true V=1.0).
    y_frac = pos0[:, 1] / L
    d_frac = 4.0 * h_mean_analytic / L
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
    kh_void = (is_sheared and kh > 0.05) or drho > 0.10

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
            band_Dx, band_Sxx, band_E = {}, {}, {}
            for label, in_blob, sign, centre in (
                ("A", in_A, +1.0, centre_A),
                ("B", in_B, -1.0, centre_B),
            ):
                x_ref = centre + np.array([sign * (v_shear / 2.0) * t, 0.0, 0.0])
                dx = min_image(pos[in_blob] - x_ref, L)[:, 0]
                w_raw = mass0[in_blob] * u[in_blob]
                w_clip = mass0[in_blob] * np.maximum(u[in_blob], 0.0)
                w_abs_sum = float(np.sum(np.abs(w_raw)))
                if w_abs_sum > 0:
                    neg_share = float(-np.sum(w_raw[w_raw < 0]) / w_abs_sum)
                    if neg_share > NEG_WEIGHT_VOID_THRESHOLD:
                        neg_weight_void = True
                Dx, Sxx = moments(w_clip, dx)
                band_Dx[label], band_Sxx[label] = Dx, Sxx
                band_E[label] = float(w_raw.sum())
            denom_s = band_Sxx["A"] + band_Sxx["B"]
            denom_e = band_E["A"] + band_E["B"]
            series[band]["A_centroid"].append(
                (band_Dx["A"] + band_Dx["B"]) / h_med_last if h_med_last > 0 else np.nan
            )
            series[band]["A_spread"].append(
                2 * (band_Sxx["A"] - band_Sxx["B"]) / denom_s if denom_s else np.nan
            )
            series[band]["A_energy"].append(
                2 * (band_E["A"] - band_E["B"]) / denom_e if denom_e else np.nan
            )

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
        sel = np.array([len(times) - 1])
    vals = values[sel]
    if not np.all(np.isfinite(vals)):
        return np.nan, int(sel.sum())
    return float(np.mean(vals)), int(sel.sum())


def group_stats(per_realisation_means):
    vals = np.array([v for v in per_realisation_means if np.isfinite(v)])
    n = len(vals)
    mean = float(np.mean(vals)) if n >= 1 else np.nan
    se = float(np.std(vals, ddof=1) / np.sqrt(n)) if n >= 2 else np.nan
    return mean, se, n


def main():
    opt = parse_options()

    shear_series = [realisation_time_series(r) for r in opt.shear_runs]
    control_series = [realisation_time_series(r) for r in opt.control_runs]

    for label, group in (("shear", shear_series), ("control", control_series)):
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
            shear_means, control_means = [], []
            for r in shear_series:
                if r["void"]:
                    continue
                avg, n_used = time_average(
                    r["times"], r["series"][band][metric], opt.window_frac
                )
                # Sign-align: v_shear < 0 mirrors the A<->B labelling.
                sign = 1.0 if r["v_shear"] >= 0 else -1.0
                shear_means.append(sign * avg)
            for r in control_series:
                if r["void"]:
                    continue
                avg, n_used = time_average(
                    r["times"], r["series"][band][metric], opt.window_frac
                )
                control_means.append(avg)

            mean_s, se_s, n_s = group_stats(shear_means)
            mean_c, se_c, n_c = group_stats(control_means)
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
            results["bands"][band][metric] = dict(
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
    sys.exit(0 if overall in ("PASS", "INCONCLUSIVE") else 1)


if __name__ == "__main__":
    main()
