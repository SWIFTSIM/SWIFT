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
Cross-run gate for ISRFMultiBinDissipation, applying the multi-bin
dissipation test specification's Sec 6 criteria. Each run's own
isrf_multibin_dissipation_check.py has already produced multibin_metrics.json;
this script only combines them.

Run roles are passed by explicit flag rather than positionally, since Sec 6's
criteria are role-specific (S0/M0/M1/M2/M3/M4/M5 each play a different part;
`--runs`/`--control` alone cannot express that unambiguously): `--s0`, `--m0`,
`--m1`, `--m2`, `--m3` (optional, interpretation aid only), `--m4`, `--m5`.
`--mode report` skips the gate and just prints every metrics file passed via
`--runs`.
"""

import argparse
import json
import os
import sys

import numpy as np

BANDS = ("PE", "LW")

# crit4's ledger bar; same derivation and value as
# isrf_multibin_dissipation_check.py's LEDGER_BAR (copied per this family's
# no-shared-module convention, not imported).
LEDGER_BAR = 1e-3


def parse_options():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--mode", choices=["gate", "report"], default="gate")
    parser.add_argument(
        "--runs", nargs="+", help="report mode: run directories to print."
    )
    parser.add_argument("--s0", help="gate mode: S0 smoke run directory.")
    parser.add_argument("--m0", help="gate mode: M0 uniform control directory.")
    parser.add_argument("--m1", help="gate mode: M1 test directory.")
    parser.add_argument(
        "--m2", help="gate mode: M2 alpha_max=0 attribution twin directory."
    )
    parser.add_argument(
        "--m3",
        help="gate mode: M3 c_hyp-pin mechanism control (interpretation aid only).",
    )
    parser.add_argument(
        "--m4", help="gate mode: M4 2x-finer-h uniform control directory."
    )
    parser.add_argument("--m5", help="gate mode: M5 2x-finer-h test directory.")
    return parser.parse_args()


def load(run_dir):
    with open(os.path.join(run_dir, "multibin_metrics.json")) as f:
        return json.load(f)


def last_third_points(series, t_max, t0=0.0):
    """Every point in the series' last third, with no filter of any kind."""
    cutoff = t0 + (2.0 / 3.0) * (t_max - t0)
    return [p for p in series if p["time"] >= cutoff]


def reduce_component(metrics, band, component, names=None):
    """Reduce one dipole component the single way every side of every
    comparison in this file uses: absolute value, over the last third of the
    series, void points dropped.

    A component recorded as None is a dipole the per-run check could not
    form, i.e. a non-finite moment. It is COUNTED, never filtered away, so
    the caller fails on it instead of quietly reducing a smaller sample.

    Parameters
    ----------
    metrics : dict
        One run's `multibin_metrics.json`.
    band : str
        "PE" or "LW".
    component : str
        Key of the component to reduce, e.g. "d_x_h".
    names : list of str, optional
        Restrict to these stars; every star by default.

    Returns
    -------
    values : list of float
    n_non_finite : int
    """
    values, n_non_finite = [], 0
    for name, per_band in metrics["per_star_series"].items():
        if names is not None and name not in names:
            continue
        series = per_band.get(band, [])
        if not series:
            continue
        t_max = max(p["time"] for p in series)
        for p in last_third_points(series, t_max):
            if p.get("void"):
                continue
            value = p.get(component)
            if value is None or not np.isfinite(value):
                n_non_finite += 1
                continue
            values.append(abs(float(value)))
    return values, n_non_finite


def compute_sigma(m0):
    """M0's own spread of |d_x|/h across its stars and last-third snapshots,
    per band: the honest within-family noise floor Sec 6 defines."""
    sigma, non_finite = {}, {}
    for band in BANDS:
        vals, n_bad = reduce_component(m0, band, "d_x_h")
        non_finite[band] = n_bad
        sigma[band] = float(np.std(vals)) if len(vals) > 1 else float("nan")
    return sigma, non_finite


def star_last(metrics, name, band):
    series = metrics["per_star_series"].get(name, {}).get(band, [])
    non_void = [p for p in series if not p.get("void")]
    return non_void[-1] if non_void else (series[-1] if series else None)


def all_points(metrics, band):
    out = []
    for name, per_band in metrics["per_star_series"].items():
        out += per_band.get(band, [])
    return out


def check_preconditions(s0, m0, sigma, sigma_non_finite):
    # Resolving power is checked PER BAND: PE and LW have genuinely
    # different physics (kappa, lambda, predicted dipole magnitude per
    # Sec 1.3), so one band failing does not make the other band's result
    # uninterpretable: it means STOP-3 fires for that band specifically
    # while the other band's verdict still stands.
    band_ok = {}
    print("=== Sec 6.1 resolving-power precondition ===")
    for band in BANDS:
        # Both sides of `med_dy >= 5*sigma` are reduced by the SAME call:
        # absolute value, last third of the series, void points dropped. The
        # left side used to keep void points and every snapshot while the
        # right side dropped void points and kept only the last third, which
        # biased the precondition towards PASS, void points being the
        # ill-conditioned ones by definition.
        d_y, n_bad_dy = reduce_component(m0, band, "d_y_h", names=["B"])
        n_bad = n_bad_dy + sigma_non_finite[band]
        if n_bad:
            print(
                f"{band}: {n_bad} non-finite dipole components inside the "
                "reduced window -> FAIL"
            )
            band_ok[band] = False
            continue
        if not d_y or not np.isfinite(sigma[band]):
            print(f"{band}: insufficient data for precondition 1 -> FAIL")
            band_ok[band] = False
            continue
        med_dy = float(np.median(d_y))
        thresh = 5.0 * sigma[band]
        this_ok = med_dy >= thresh
        band_ok[band] = this_ok
        print(
            f"{band}: median|d_y/h| at B (M0) = {med_dy:.5f}  "
            f"5*sigma = {thresh:.5f}  sigma={sigma[band]:.5f}  "
            f"-> {'PASS' if this_ok else 'FAIL (insufficient resolving power) -- STOP-3 for this band'}"
        )
    bin_span_ok = s0["bin_span"] >= 2.0
    phases_ok = len(s0.get("drho_by_phase", {})) >= 2
    print(
        f"S0 bin span = {s0['bin_span']:.3f}  (need >= 2)  "
        f"phases populated = {list(s0.get('drho_by_phase', {}).keys())}  "
        f"-> {'PASS' if (bin_span_ok and phases_ok) else 'FAIL (no bin split)'}"
    )
    if not (bin_span_ok and phases_ok):
        band_ok = {band: False for band in BANDS}

    v_ratio = s0.get("v_over_c_hyp_realized")
    if v_ratio is not None:
        in_window = 0.18 <= v_ratio <= 0.32
        print(
            f"S0 realized v/c_hyp = {v_ratio:.4f}  "
            f"(window [0.18, 0.32])  -> {'ok' if in_window else 'STOP-1b: re-derive star velocity and re-run S0'}"
        )
    return band_ok


def mode_gate(opt):
    if not all([opt.s0, opt.m0, opt.m1, opt.m2, opt.m4, opt.m5]):
        raise RuntimeError("gate mode requires --s0 --m0 --m1 --m2 --m4 --m5")

    s0, m0, m1, m2, m4, m5 = (
        load(opt.s0),
        load(opt.m0),
        load(opt.m1),
        load(opt.m2),
        load(opt.m4),
        load(opt.m5),
    )
    m3 = load(opt.m3) if opt.m3 else None

    excluded = []
    for tag, m in (("M0", m0), ("M1", m1), ("M2", m2), ("M4", m4), ("M5", m5)):
        if m.get("invalid_bin_dt") or m.get("void_drho") or m.get("void_dx_gas"):
            excluded.append(tag)
    if excluded:
        print(f"Excluded from the gate (VOID/INVALID): {excluded}")

    sigma, sigma_non_finite = compute_sigma(m0)
    print(f"sigma (M0 |d_x|/h spread): {sigma}")

    band_ok = check_preconditions(s0, m0, sigma, sigma_non_finite)
    if not any(band_ok.values()):
        print(
            "\nOverall gate: INVALID -- resolving-power precondition failed for every band. No verdict."
        )
        sys.exit(1)

    if not m1.get("ab_matched", True):
        print(
            f"\nSTOP-2: M1 rho_match = {m1.get('rho_match')} exceeds 0.15 -- "
            "star A and star B are not a matched pair. crit1 (the within-run "
            "matched-pair criterion) is reported below but is not a valid "
            "primary verdict; fall back to interpreting M1 against the M0 "
            "control per Sec 3.4's own designated fallback for this case."
        )

    overall_ok = True
    print("\n=== Sec 6.2 verdict ===")
    for band in BANDS:
        print(f"-- {band} --")
        if not band_ok[band]:
            print(
                f"  INVALID -- resolving-power precondition failed for {band} (STOP-3). No verdict for this band."
            )
            continue
        A1 = star_last(m1, "A", band)
        B1 = star_last(m1, "B", band)
        if A1 is None or B1 is None:
            print("  missing A/B data in M1 -> FAIL")
            overall_ok = False
            continue
        # np.maximum, not the builtin: `max(a, nan)` returns `a`, so a
        # non-finite sibling dipole would set the limit to 3*sigma and the
        # criterion would pass on evidence it never saw.
        a_dx = A1["d_x_h"]
        b_dx = B1["d_x_h"]
        if a_dx is None or b_dx is None:
            print("  crit1: a non-finite dipole at A or B -> FAIL")
            overall_ok = False
            continue
        limit1 = float(np.maximum(abs(b_dx), 3.0 * sigma[band]))
        crit1 = bool(np.isfinite(limit1) and abs(a_dx) <= limit1)
        print(
            f"  crit1 (matched pair): |d_x/h|@A={abs(a_dx):.5f}  "
            f"limit=max(|d_x/h|@B={abs(b_dx):.5f}, 3*sigma={3*sigma[band]:.5f})="
            f"{limit1:.5f}  -> {'PASS' if crit1 else 'FAIL'}"
        )
        overall_ok &= crit1

        # VOID points are excluded from every gated metric (Sec 3.4), not
        # just crit1's matched-pair criterion; a VOID star's own >0.1h is
        # still printed (not silently dropped) since it can be a real
        # finding even though it does not gate.
        crit2 = True
        for tag, m in (("M0", m0), ("M1", m1), ("M2", m2), ("M4", m4), ("M5", m5)):
            for p in all_points(m, band):
                if p.get("d_total_h") is None:
                    if not p.get("void"):
                        crit2 = False
                        print(
                            f"  crit2 VIOLATED in {tag}: non-finite |d|/h at "
                            f"t={p['time']}"
                        )
                    continue
                if p["d_total_h"] <= 0.1:
                    continue
                if p.get("void"):
                    print(
                        f"  crit2 note (VOID, not gated) in {tag}: |d|/h = {p['d_total_h']:.5f} > 0.1 at t={p['time']}"
                    )
                    continue
                crit2 = False
                print(
                    f"  crit2 VIOLATED in {tag}: |d|/h = {p['d_total_h']:.5f} > 0.1 at t={p['time']}"
                )
        print(
            f"  crit2 (|d|/h <= 0.1 everywhere, non-void): -> {'PASS' if crit2 else 'FAIL'}"
        )
        overall_ok &= crit2

        A5 = star_last(m5, "A", band)
        crit3 = True
        if (
            A1
            and A5
            and A1.get("d_total_h") is not None
            and A5.get("d_total_h") is not None
        ):
            growth = A5["d_total_h"] - A1["d_total_h"]
            crit3 = bool(
                np.isfinite(growth)
                and np.isfinite(sigma[band])
                and growth <= 3.0 * sigma[band]
            )
            print(
                f"  crit3 (no growth with resolution): |d|/h@A M1={A1['d_total_h']:.5f} "
                f"M5={A5['d_total_h']:.5f}  growth={growth:.5f}  3*sigma={3*sigma[band]:.5f}  "
                f"-> {'PASS' if crit3 else 'FAIL'}"
            )
        else:
            print("  crit3: missing data -> FAIL")
            crit3 = False
        overall_ok &= crit3

        # crit4: the closed ledger residual |E + Abs - Inj| / |Inj|, the same
        # reduction ISRFInjectionConservation/isrf_ledger_check.py gates.
        # This is the only sound reference for Sec 6.2's non-conservation
        # question (see REVIEW_multibin_crit4_and_opticallythin_window):
        # `fractional_change`, compared previously, is the field's GROWTH
        # under continuous injection, not a conservation residual, and M0/M1
        # differ in density structure so their growth curves differ for
        # reasons unrelated to conservation. The residual is absolute and
        # zero under exact conservation for any IC, so M0 (different IC) is
        # a legitimate "do seams leak" reference; M2 (M1's alpha_max=0 twin
        # on the IDENTICAL IC) is reported alongside as the attribution half.
        #
        # The two accumulators this needs are identically zero unless swift
        # was built with --enable-debugging-checks (feedback_common.c:1248);
        # run.sh detects this from the binary's own --version banner and
        # passes it through as ledger_valid. A production-build campaign
        # SKIPs crit4 rather than gating on all-zero data.
        ledger_ok = all(m.get("ledger_valid") for m in (m0, m1, m2) if m is not None)
        if not ledger_ok:
            print(
                "  crit4 (ledger residual, M1 vs M0): SKIPPED, not a "
                "debugging-build campaign (ledger_valid=False on at least "
                "one of M0/M1/M2). Not gated; does not fail the run."
            )
        else:
            R0 = m0["m4_ledger"][band]["max_R"]
            R1 = m1["m4_ledger"][band]["max_R"]
            R2 = m2["m4_ledger"][band]["max_R"]
            limit4 = float(np.maximum(R0, LEDGER_BAR)) if np.isfinite(R0) else np.nan
            crit4 = bool(np.isfinite(R1) and np.isfinite(limit4) and R1 <= limit4)
            print(
                f"  crit4 (ledger residual, M1 vs M0): max_R M0={R0:.3e}  "
                f"M1={R1:.3e}  limit=max(M0, {LEDGER_BAR:g})={limit4:.3e}  "
                f"M2 (attribution)={R2:.3e}  -> {'PASS' if crit4 else 'FAIL'}"
            )
            overall_ok &= crit4

        # FAIL's attribution half: crit1 fails AND M2's same quantity is void of signal.
        if not crit1:
            A2 = star_last(m2, "A", band)
            if A2 is not None:
                attributable = abs(A2["d_x_h"]) < 3.0 * sigma[band]
                print(
                    f"  interpretation: M2 (alpha_max=0) |d_x/h|@A = {abs(A2['d_x_h']):.5f}  "
                    f"3*sigma={3*sigma[band]:.5f}  -> "
                    f"{'attributable to dissipation' if attributable else 'NOT attributable to dissipation (interface artifact) -- STOP-5'}"
                )

    if m3 is not None:
        print(
            "\n=== M3 mechanism control (interpretation aid, two-factor, do not over-read) ==="
        )
        for band in BANDS:
            A3 = star_last(m3, "A", band)
            if A3 is not None:
                collapsed = abs(A3["d_x_h"]) < 3.0 * sigma[band]
                print(
                    f"  {band}: |d_x/h|@A under c_hyp pin = {abs(A3['d_x_h']):.5f}  "
                    f"3*sigma={3*sigma[band]:.5f}  -> "
                    f"{'consistent with D-step mechanism' if collapsed else 'survives pin: points away from min(c_hyp)'}"
                )

    print(f"\nOverall gate: {'PASS' if overall_ok else 'FAIL'}")
    if not overall_ok:
        sys.exit(1)


def mode_report(opt):
    if not opt.runs:
        raise RuntimeError("report mode requires --runs")
    for run_dir in opt.runs:
        m = load(run_dir)
        print(f"=== {run_dir} ===")
        print(
            json.dumps(
                {
                    k: m[k]
                    for k in (
                        "variant",
                        "bin_delta",
                        "status",
                        "bin_span",
                        "v_over_c_hyp_realized",
                        "rho_match",
                        "ab_matched",
                    )
                    if k in m
                },
                indent=2,
            )
        )


def main():
    opt = parse_options()
    if opt.mode == "gate":
        mode_gate(opt)
    else:
        mode_report(opt)


if __name__ == "__main__":
    main()
