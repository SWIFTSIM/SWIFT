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
Cross-run comparisons for ISRFPropagationSpeedSweep: M-C1/M-C2 (Leg P), M-C3
(Leg N), M-C4 (Legs R and I, report only). Each run directory must already
contain its own sweep_metrics.json (isrf_propagation_speed_sweep_check.py)
and its last snapshot.
"""

import argparse
import glob
import json
import os
import sys

import h5py
import numpy as np


def parse_options():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--mode", choices=["P", "N", "RI"], required=True,
        help="Which cross-run gate to compute.",
    )
    parser.add_argument(
        "--runs", nargs="+", required=True,
        help="Run directories, in the order given by the run table (Leg P: "
        "P1..P7; Leg N: N1..N8; Leg RI: R1..R5 then I1..I4).",
    )
    parser.add_argument(
        "--ref-index", type=int, default=3,
        help="Leg P only: 0-based index of the reference run in --runs "
        "(default 3 = P4, mid-sweep).",
    )
    parser.add_argument("--tol-d", type=float, default=1e-3, help="M-C1 D tolerance.")
    parser.add_argument(
        "--tol-retardation", type=float, default=0.10, help="M-C2 relative tolerance."
    )
    parser.add_argument(
        "--bracket-factor", type=float, default=1.25,
        help="M-C3 bracket window (0.8x .. 1.25x nu_max).",
    )
    return parser.parse_args()


def load_metrics(run_dir):
    with open(os.path.join(run_dir, "sweep_metrics.json")) as f:
        return json.load(f)


def last_snapshot(run_dir):
    files = sorted(glob.glob(os.path.join(run_dir, "snap", "snapshot_*.hdf5")))
    if not files:
        raise RuntimeError(f"No snapshots found in {run_dir}/snap")
    with h5py.File(files[-1], "r") as f:
        gas = f["/PartType0"]
        ids = gas["ParticleIDs"][:]
        u_fuv = gas["FUVSpecificEnergies"][:].astype(np.float64)
        u_lw = gas["LWSpecificEnergies"][:].astype(np.float64)
    return ids, u_fuv, u_lw


def mode_P(opt):
    metrics = [load_metrics(d) for d in opt.runs]

    # M-P6 validity precondition: nu_eff spread and step-count agreement.
    nu_effs = np.array([m["nu_eff"] for m in metrics])
    n_steps = np.array([m["n_steps"] for m in metrics])
    nu_spread = (nu_effs.max() - nu_effs.min()) / np.mean(nu_effs) if np.mean(nu_effs) else np.inf
    steps_agree = len(set(n_steps.tolist())) == 1
    valid = (nu_spread <= 0.01) and steps_agree

    print("=== Leg P: M-P6 validity precondition ===")
    print(f"nu_eff per run: {nu_effs}")
    print(f"step counts per run: {n_steps}")
    print(f"nu_eff spread: {nu_spread:.4%}  (limit 1%)")
    print(f"step counts agree: {steps_agree}")
    print(f"=> {'VALID' if valid else 'INVALID'}")

    if not valid:
        print("\nM-C1 and M-C2: INVALID (not FAILED) -- the dt-quantization "
              "precondition is not met. Re-run with time_end adjusted so the "
              "bin grid lands on the intended dt; do not interpret D below.")

    drift_void = np.array([m["drift_void"] for m in metrics])

    # M-C1: self-similarity D against the reference run.
    ref_ids, ref_u_fuv, ref_u_lw = last_snapshot(opt.runs[opt.ref_index])
    ref_sort = np.argsort(ref_ids)
    print("\n=== M-C1: self-similarity (Leg P) ===")
    all_pass = True
    for i, d in enumerate(opt.runs):
        ids, u_fuv, u_lw = last_snapshot(d)
        order = np.argsort(ids)
        if not np.array_equal(ids[order], ref_ids[ref_sort]):
            print(f"{d}: ParticleIDs do not match the reference run -- skipped.")
            continue
        d_fuv = np.max(np.abs(u_fuv[order] - ref_u_fuv[ref_sort])) / max(
            np.max(np.abs(ref_u_fuv)), 1e-300
        )
        d_lw = np.max(np.abs(u_lw[order] - ref_u_lw[ref_sort])) / max(
            np.max(np.abs(ref_u_lw)), 1e-300
        )
        excluded = bool(drift_void[i])
        status = "EXCLUDED (drift > 0.1h)" if excluded else (
            "PASS" if valid and max(d_fuv, d_lw) <= opt.tol_d else "FAIL"
        )
        if not excluded and valid and max(d_fuv, d_lw) > opt.tol_d:
            all_pass = False
        print(f"{d}: D_FUV={d_fuv:.3e}  D_LW={d_lw:.3e}  -> {status}")
    print(f"M-C1 overall: {'PASS' if (valid and all_pass) else ('INVALID' if not valid else 'FAIL')}")

    # M-C2: retardation collapse at eps=0.01.
    print("\n=== M-C2: retardation collapse (Leg P, eps=0.01) ===")
    norms_fuv, norms_lw = [], []
    for i, (d, m) in enumerate(zip(opt.runs, metrics)):
        r_fuv = m["front"]["FUV"]["0.01"]["r_edge_norm"] if "0.01" in m["front"]["FUV"] else m["front"]["FUV"][0.01]["r_edge_norm"]
        r_lw = m["front"]["LW"]["0.01"]["r_edge_norm"] if "0.01" in m["front"]["LW"] else m["front"]["LW"][0.01]["r_edge_norm"]
        norms_fuv.append(r_fuv)
        norms_lw.append(r_lw)
        void_str = " (drift-void, still reported)" if drift_void[i] else ""
        print(f"{d}: r_edge/(c_hyp*t) FUV={r_fuv:.4f}  LW={r_lw:.4f}{void_str}")
    kept_fuv = [v for v, void in zip(norms_fuv, drift_void) if not void]
    kept_lw = [v for v, void in zip(norms_lw, drift_void) if not void]
    spread_fuv = (max(kept_fuv) - min(kept_fuv)) / np.mean(kept_fuv) if kept_fuv else np.inf
    spread_lw = (max(kept_lw) - min(kept_lw)) / np.mean(kept_lw) if kept_lw else np.inf
    c2_pass = valid and spread_fuv <= opt.tol_retardation and spread_lw <= opt.tol_retardation
    print(f"Spread (drift-valid runs only): FUV={spread_fuv:.4%}  LW={spread_lw:.4%}  (limit {opt.tol_retardation:.0%})")
    print(f"M-C2: {'PASS' if c2_pass else ('INVALID' if not valid else 'FAIL')}")


def mode_N(opt):
    metrics = [load_metrics(d) for d in opt.runs]
    print("=== Leg N: M-C3 bound location ===")
    stable = [m["stability_ok"] for m in metrics]
    nu_eff = [m["nu_eff"] for m in metrics]
    nu_max = [m["nu_max"] for m in metrics]
    for d, s, n, nm in zip(opt.runs, stable, nu_eff, nu_max):
        print(f"{d}: nu_eff={n:.4f}  nu_max(predicted)={nm:.4f}  "
              f"stable={s}")
    # Group by alpha (nu_max value) since N7/N8 use alpha=0.
    by_alpha = {}
    for d, s, n, nm in zip(opt.runs, stable, nu_eff, nu_max):
        by_alpha.setdefault(round(nm, 4), []).append((d, s, n))
    overall_pass = True
    for nm, entries in by_alpha.items():
        stable_nu = [n for (_, s, n) in entries if s]
        unstable_nu = [n for (_, s, n) in entries if not s]
        largest_stable = max(stable_nu) if stable_nu else None
        smallest_unstable = min(unstable_nu) if unstable_nu else None
        print(f"\n-- Group nu_max={nm:.4f} --")
        print(f"largest stable nu_eff: {largest_stable}")
        print(f"smallest unstable nu_eff: {smallest_unstable}")
        bracket_ok = True
        for (d, s, n) in entries:
            if not s and n < 0.8 * nm:
                print(f"  FAIL: {d} at nu_eff={n} < 0.8*nu_max={0.8*nm:.4f} is unstable "
                      f"-- bound is optimistic.")
                bracket_ok = False
            if s and n > opt.bracket_factor * nm:
                print(f"  FAIL: {d} at nu_eff={n} > {opt.bracket_factor}*nu_max="
                      f"{opt.bracket_factor*nm:.4f} is stable -- bound over-restricts.")
                bracket_ok = False
        if largest_stable is not None and smallest_unstable is not None:
            brackets = largest_stable < nm < smallest_unstable or \
                (largest_stable <= nm and smallest_unstable >= nm)
            print(f"  Bracket [{largest_stable}, {smallest_unstable}] contains nu_max={nm:.4f}: {brackets}")
        overall_pass &= bracket_ok
    print(f"\nM-C3 overall: {'PASS' if overall_pass else 'FAIL'}")
    if not overall_pass:
        sys.exit(1)


def mode_RI(opt):
    metrics = [load_metrics(d) for d in opt.runs]
    print("=== M-C4: Legs R and I (report only, no gate) ===")
    SPEED_OF_LIGHT_KM_S = 2.99792458e5
    for d, m in zip(opt.runs, metrics):
        c_hyp = m["c_hyp"]
        print(f"{d}: c_hyp={c_hyp:.4f} km/s  c/c_hyp={SPEED_OF_LIGHT_KM_S / c_hyp if c_hyp > 0 else float('nan'):.3e}")


def main():
    opt = parse_options()
    if opt.mode == "P":
        mode_P(opt)
    elif opt.mode == "N":
        mode_N(opt)
    else:
        mode_RI(opt)


if __name__ == "__main__":
    main()
