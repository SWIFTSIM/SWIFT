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
Relative-velocity lag of the ISRF field behind a star moving through static
gas (theory chapter, `sec:pe-lagrangian` of `02_fuv_isrf.tex`).

Measured lag, per band: the mass-weighted first moment of u over the whole
periodic box, relative to the star, along its direction of motion,
  lag = -sum_i m_i u_i xi_i / sum_i m_i u_i,  xi_i = minimum image of x_i - x_star,
with no radius window and no u > 0 filter. Analytic prediction for a
uniform c_hyp and no dissipation, from the steady moment balance of the
continuum equation: lag = v_rel tau a / (1 - exp(-a)), tau = lambda/c_hyp,
a = dt/tau, lambda = 1/(kappa rho) from the code's dust opacity.

That prediction already assumes a *steady* field: sum_i m_i xi_i (div F)_i
is not zero in general (it equals -integral rho F over the box), only once
the net flux itself has relaxed to zero, which is what "steady" below
checks for. The discrete estimator (`{band}SpecificFluxDivergences`,
`(1/rho) div(rho F)` per particle, finalized in the density ghost)
conserves the *zeroth* moment sum_i m_i (div F)_i = 0 exactly, because the
pairwise contribution to one particle is minus the contribution to its
neighbour (`sec:pe-operators`, `radiation_propagation_iact.h`:
`div_F_i += mj*Phi_ij`, `div_F_j += -mi*Phi_ij`). Its *first* moment does
not cancel the same way (a pair's two xi differ), so the discrete scheme
carries its own first-moment residual, a net-flux residual plausibly
sourced by the propagation operator's non-antisymmetric own-derivative
gradient term, that the steady-state prediction above omits:
  T2 = tau * sum_i m_i xi_i (div F)_i / sum_i m_i u_i,
same units as `lag` (tau times a specific-energy-rate first moment, divided
by a specific-energy zeroth moment, is a length). T2 IS NOT GATED and never
was a reference: it is measured from the run itself, from a field the
propagation operator under test writes, so `pred + T2` absorbs whatever
net-flux residual is actually present and cannot discriminate it. It is
printed for comparison, together with `lag/(pred + T2) - 1`.

The gate uses `lag/pred - 1`, whose reference is independent of the u field:
`tau = lambda/c_hyp` with lambda from the dust constants re-derived here,
`c_hyp` from `used_parameters.yml` and `v_rel` from the ICs.

Bar. The terms of `|lag/pred - 1|` that CAN be derived in advance are the
minimum-image wrap of the exponential profile, `exp(-L/(2 lambda))`, which is
below 1e-6 at the required `L >= 30 lambda`; the float32 storage of u inside
the first moment; and the glass's own density scatter entering `lambda_med`,
a few 1e-3. They are printed per snapshot and sum well below REL_BAR. What
CANNOT be derived is the discrete estimator's own first-moment residual: this
fixture runs at `lambda_med/h_med` of about 1 to 1.6, so the absorption
length is not resolved by the kernel and no a-priori coefficient exists for
that term. T2 is a measurement of it, not a derivation, which is exactly why
it may not set the bar. REL_BAR = 0.02 is therefore a DECLARED bar sitting
above the derivable floor, and this limit is a live open item, not a
derivation.

--gate (alpha_max = alpha_floor = 0 required): per band, evaluate
`|lag/pred - 1| <= 0.02` on every snapshot the field is STEADY,
where steady is decided *before* looking at the ratio and requires BOTH:
(1) the relative change of sum_i m_i u_i between that snapshot and the
previous one is < 1e-3 (the zeroth moment has stopped changing), and
(2) the relative change of `lag` itself (the first moment used by this
check) between that snapshot and the previous one is < LAG_REL_TOL, which
is half of REL_BAR: a snapshot whose own last step still moved the gated
numerator by more than half the bar is not converged.
Criterion (1) alone is not enough: the zeroth moment can settle while
`lag` itself is still oscillating step to step, or while the field's
*shape* (not just its total energy) is still settling.
INCONCLUSIVE exits 1, like FAIL: a gate that could not be evaluated is not
a pass.
PASS if every doubly-steady snapshot of both bands passes; INCONCLUSIVE if
neither band ever FAILs but a band has no doubly-steady snapshot in the
evaluated window (the gate cannot be evaluated, not a failure; distinct
from INVALID below); FAIL on any non-finite value or a doubly-steady
snapshot outside 2%. A band with no qualifying snapshot is always reported
(with a geometric-decay estimate of the extra run length needed, when the
trend allows one), even when the other band's own violation already made
the run FAIL. The old analytic-only ratio `lag/pred - 1` (missing T2) is
still printed, informational only, together with the T2-corrected ratio.

--report-only (dissipation on): print the displacement (lag - pred - T2)/h
next to the 0.05 h reference, exit 0. Any non-finite value: FAIL. The box
must satisfy L >= 30 lambda in both bands for the gate; a smaller box is
reported as INVALID (a geometry failure, not the INCONCLUSIVE above).
"""

import argparse
import glob
import json
import os
import sys

import h5py
import numpy as np
import yaml

# src/feedback/GEAR/radiation.h
SIGMA_D_CGS = {"PE": 9e-22, "LW": 1.5e-21}
MU_H = 1.4
M_H_CGS = 1.6726219e-24
# Grackle's SolarMetalFractionByMass default
GRACKLE_SOLAR_Z = 0.01295
REL_BAR = 0.02
H_REFERENCE = 0.05
MIN_BOX_OVER_LAMBDA = 30.0
STEADY_REL_TOL = 1e-3
# A snapshot whose own last step moved `lag` by more than half the bar cannot
# be called converged: that step's unconverged remainder alone would use up
# half the budget the gated ratio is judged against. Derived from REL_BAR, not
# from any observed value.
LAG_REL_TOL = 0.5 * REL_BAR
# Float32 storage of u inside the mass-weighted first moment.
FLOAT32_EPS = 1.19209e-07
# How many trailing snapshots to evaluate the gate/report over (needs one
# extra leading snapshot per evaluated one, to test its own steadiness).
N_EVAL = 3


def parse_options():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--run", type=str, required=True, help="Run directory.")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--gate", action="store_true")
    mode.add_argument("--report-only", action="store_true")
    return parser.parse_args()


def modal_dt(run):
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


def lag_of(path, v_hat):
    with h5py.File(path, "r") as f:
        units = f["/Units"].attrs
        ul = float(np.asarray(units["Unit length in cgs (U_L)"]).flat[0])
        um = float(np.asarray(units["Unit mass in cgs (U_M)"]).flat[0])
        L = float(np.asarray(f["/Header"].attrs["BoxSize"]).flat[0])
        t = float(np.asarray(f["/Header"].attrs["Time"]).flat[0])
        gas = f["/PartType0"]
        pos = gas["Coordinates"][:].astype(np.float64)
        m = gas["Masses"][:].astype(np.float64)
        rho = gas["Densities"][:].astype(np.float64)
        h = gas["SmoothingLengths"][:].astype(np.float64)
        Z = gas["MetalMassFractions"][:, -1].astype(np.float64)
        u = {b: gas[f"{b}SpecificEnergies"][:].astype(np.float64) for b in SIGMA_D_CGS}
        div_f = {
            b: gas[f"{b}SpecificFluxDivergences"][:].astype(np.float64)
            for b in SIGMA_D_CGS
        }
        star = f["/PartType4/Coordinates"][0, :].astype(np.float64)
    xi = pos - star
    xi -= L * np.round(xi / L)
    xi_par = xi @ v_hat
    out = dict(time=t, boxsize=L, h_med=float(np.median(h)), bands={})
    for band, sigma in SIGMA_D_CGS.items():
        kappa = sigma * (Z / GRACKLE_SOLAR_Z) / (MU_H * M_H_CGS)
        lam = 1.0 / (kappa * rho * um / ul**3) / ul
        finite = bool(np.all(np.isfinite(u[band]))) and bool(
            np.all(np.isfinite(div_f[band]))
        )
        M0 = float(np.sum(m * u[band]))
        lag = (
            float(-np.sum(m * u[band] * xi_par) / M0) if finite and M0 != 0 else np.nan
        )
        # First moment of the discrete divergence accumulator, along xi_par;
        # tau (band- and closure-dependent) is applied by the caller, once
        # c_hyp is known.
        first_moment_div_f = (
            float(np.sum(m * xi_par * div_f[band]) / M0)
            if finite and M0 != 0
            else np.nan
        )
        lam_med = float(np.median(lam))
        out["bands"][band] = dict(
            lag=lag,
            lambda_med=lam_med,
            n_gas=int(m.size),
            # lambda goes as 1/rho, so the glass's own density scatter sets
            # how far the median lambda the prediction uses sits from the
            # lambda the first moment actually samples.
            lambda_scatter=(
                float(np.median(np.abs(lam - lam_med)) / lam_med)
                if lam_med > 0
                else np.inf
            ),
            finite=finite,
            M0=M0,
            first_moment_div_f=first_moment_div_f,
            negative_mass_share=(
                float(np.sum(m * np.minimum(u[band], 0)) / M0) if M0 != 0 else np.nan
            ),
        )
    return out


def _extrapolate_series(times, values, tol):
    """Geometric-decay extrapolation: from the last two finite, decreasing
    values below tol's reach, how many more equally-spaced points until the
    series drops below tol. None if fewer than 2 usable points, not
    decreasing, or already converged."""
    pts = [(t, v) for t, v in zip(times, values) if np.isfinite(v)]
    if len(pts) < 2:
        return None
    (t0, v0), (t1, v1) = pts[-2], pts[-1]
    if not (v1 < v0) or v1 <= 0 or v1 <= tol:
        return None
    rate = v1 / v0
    dt_snap = t1 - t0
    if dt_snap <= 0:
        return None
    n_more = int(np.ceil(np.log(tol / v1) / np.log(rate)))
    return n_more, dt_snap, t1


def estimate_more_snapshots(history):
    """history: list of (time, rel_M0, rel_lag). Tries both criteria's own
    trend and returns the one needing MORE additional snapshots (the
    binding constraint, since both must hold together), as
    (label, n_more, dt_snap, t_last, t_extra); None if neither extrapolates."""
    times = [h[0] for h in history]
    candidates = []
    for label, tol, col in (("rel_M0", STEADY_REL_TOL, 1), ("rel_lag", LAG_REL_TOL, 2)):
        values = [h[col] for h in history]
        r = _extrapolate_series(times, values, tol)
        if r is not None:
            n_more, dt_snap, t_last = r
            candidates.append((label, n_more, dt_snap, t_last))
    if not candidates:
        return None
    label, n_more, dt_snap, t_last = max(candidates, key=lambda c: c[1])
    return label, n_more, dt_snap, t_last, n_more * dt_snap


def main():
    opt = parse_options()
    params = yaml.safe_load(open(os.path.join(opt.run, "used_parameters.yml")))
    fb = params["GEARFeedback"]
    c_hyp = float(fb["ISRF_c_hyp_pin_for_debugging"])
    if c_hyp <= 0.0:
        print("FAIL: this check needs ISRF_c_hyp_pin_for_debugging > 0 (uniform tau).")
        sys.exit(1)
    alpha_max = float(fb["ISRF_dissipation_alpha_max"])
    alpha_floor = float(fb["ISRF_dissipation_alpha_floor"])
    if opt.gate and (alpha_max != 0.0 or alpha_floor != 0.0):
        print(
            "FAIL: --gate requires ISRF_dissipation_alpha_max = "
            f"ISRF_dissipation_alpha_floor = 0 (got {alpha_max}, {alpha_floor}); "
            "the 2% T2-corrected bar is only adjudicated at alpha 0, use "
            "--report-only for a dissipative run."
        )
        sys.exit(1)
    with h5py.File(
        os.path.join(opt.run, "ICs_isrf_galilean_invariance.hdf5"), "r"
    ) as f:
        v_star = f["/PartType4/Velocities"][0, :].astype(np.float64)
        v_gas = f["/PartType0/Velocities"][:].astype(np.float64).mean(axis=0)
    v_rel_vec = v_star - v_gas
    v_rel = float(np.linalg.norm(v_rel_vec))
    v_hat = v_rel_vec / v_rel
    dt = modal_dt(opt.run)
    files = sorted(glob.glob(os.path.join(opt.run, "snap", "snapshot_*.hdf5")))
    if not files:
        print(f"FAIL: no snapshots found in {opt.run}/snap/snapshot_*.hdf5")
        sys.exit(1)
    print(
        f"{opt.run}: v_rel={v_rel:.4g} km/s = {v_rel / c_hyp:.3g} c_hyp, c_hyp={c_hyp} km/s, "
        f"dt={dt:.4e}, alpha_max={alpha_max}, alpha_floor={alpha_floor}"
    )

    # One extra leading snapshot so the first evaluated one can test its own
    # steadiness against its predecessor.
    window = files[-(N_EVAL + 1) :] if len(files) > N_EVAL else files
    results = [lag_of(path, v_hat) for path in window]

    verdict = "PASS"
    any_steady = {b: False for b in SIGMA_D_CGS}
    rel_history = {b: [] for b in SIGMA_D_CGS}  # (time, rel_M0, rel_lag) per band
    summary = dict(
        run=opt.run,
        v_rel_kms=v_rel,
        c_hyp_kms=c_hyp,
        dt=dt,
        steady_rel_tol=STEADY_REL_TOL,
        lag_rel_tol=LAG_REL_TOL,
        snapshots={},
    )
    prev = None
    for idx, (path, res) in enumerate(zip(window, results)):
        eval_this = idx >= len(window) - N_EVAL
        for band, r in res["bands"].items():
            tau = r["lambda_med"] / c_hyp
            a = dt / tau
            pred = v_rel * tau * a / (-np.expm1(-a))
            T2 = tau * r["first_moment_div_f"]
            pred_corr = pred + T2
            r.update(
                tau=tau,
                a=a,
                prediction=pred,
                T2=T2,
                prediction_corrected=pred_corr,
                t_over_tau=res["time"] / tau,
                box_over_lambda=res["boxsize"] / r["lambda_med"],
            )
            ratio_old = r["lag"] / pred
            ratio_new = r["lag"] / pred_corr if pred_corr != 0 else np.nan
            disp_h = (r["lag"] - pred_corr) / res["h_med"]
            rel_M0 = rel_lag = np.nan
            steady_M0 = steady_lag = False
            if (
                prev is not None
                and np.isfinite(prev[band]["M0"])
                and prev[band]["M0"] != 0
            ):
                rel_M0 = abs(r["M0"] - prev[band]["M0"]) / abs(prev[band]["M0"])
                steady_M0 = rel_M0 < STEADY_REL_TOL
            if (
                prev is not None
                and np.isfinite(prev[band]["lag"])
                and prev[band]["lag"] != 0
                and np.isfinite(r["lag"])
            ):
                rel_lag = abs(r["lag"] - prev[band]["lag"]) / abs(prev[band]["lag"])
                steady_lag = rel_lag < LAG_REL_TOL
            steady = steady_M0 and steady_lag
            r.update(
                ratio_old=ratio_old,
                ratio_new=ratio_new,
                displacement_h=disp_h,
                rel_M0=rel_M0,
                rel_lag=rel_lag,
                steady=bool(steady),
            )
            rel_history[band].append((res["time"], rel_M0, rel_lag))
            if eval_this:
                # Derivable terms of the bar, per snapshot. The unresolved
                # kernel-scale term is NOT among them; see the docstring.
                wrap = float(np.exp(-0.5 * res["boxsize"] / r["lambda_med"]))
                storage = FLOAT32_EPS * np.sqrt(float(r["n_gas"]))
                floor = wrap + storage + r["lambda_scatter"]
                print(
                    f"  t={res['time']:.3e} {band}: derivable bar floor="
                    f"{floor:.2e} (wrap {wrap:.1e} + float32 {storage:.1e} + "
                    f"lambda scatter {r['lambda_scatter']:.1e}) vs bar {REL_BAR:g}"
                )
                print(
                    f"  t={res['time']:.3e} {band}: lag={r['lag']:.4e} pred={pred:.4e} "
                    f"T2={T2:.4e} pred+T2={pred_corr:.4e} old-ratio-1={ratio_old - 1:+.4f} "
                    f"new-ratio-1={ratio_new - 1:+.4f} disp={disp_h:+.4f} h a={a:.4f} "
                    f"t/tau={r['t_over_tau']:.1f} L/lambda={r['box_over_lambda']:.1f} "
                    f"rel_M0={rel_M0:.4f} rel_lag={rel_lag:.4f} steady={r['steady']} "
                    f"neg-share={r['negative_mass_share']:.3f}"
                )
                if not r["finite"] or not np.isfinite(ratio_old):
                    verdict = "FAIL"
                elif opt.gate:
                    if r["box_over_lambda"] < MIN_BOX_OVER_LAMBDA and verdict != "FAIL":
                        verdict = "INVALID"
                    elif r["steady"]:
                        any_steady[band] = True
                        if abs(ratio_old - 1.0) > REL_BAR:
                            verdict = "FAIL"
        prev = {band: r for band, r in res["bands"].items()}
        summary["snapshots"][os.path.basename(path)] = res
    if opt.gate:
        # Reported for any band with zero qualifying snapshots, regardless of
        # whether the other band already made the run FAIL: a per-band gap in
        # coverage is worth knowing about either way.
        for band, ok in any_steady.items():
            if ok:
                continue
            print(
                f"  {band}: no snapshot both sum-m-u-steady (< {STEADY_REL_TOL:g}) "
                f"and lag-steady (< {LAG_REL_TOL:g}) in the evaluated window."
            )
            more = estimate_more_snapshots(rel_history[band])
            if more is None:
                print(
                    "    Cannot extrapolate a run length: neither rel_M0 nor "
                    "rel_lag decays monotonically over the evaluated window; "
                    "rerun with more snapshots past the current time_end instead "
                    "of trusting an extrapolation."
                )
            else:
                which, n_more, dt_snap, t_last, t_extra = more
                print(
                    f"    Extrapolated from {which}'s geometric decay: about "
                    f"{n_more} more snapshot interval(s) of {dt_snap:.3e} each "
                    f"(+{t_extra:.3e} code time, new time_end >= "
                    f"{t_last + t_extra:.3e})."
                )
        if verdict == "PASS" and not all(any_steady.values()):
            verdict = "INCONCLUSIVE"
    if opt.report_only and verdict == "PASS":
        verdict = "REPORT"
        print(
            f"  report only: displacement against the {H_REFERENCE} h reference "
            "(T2-corrected prediction)"
        )
    summary["verdict"] = verdict
    with open(os.path.join(opt.run, "lag_metrics.json"), "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print(f"{opt.run}: lag check {verdict}")
    sys.exit(0 if verdict in ("PASS", "REPORT") else 1)


if __name__ == "__main__":
    main()
