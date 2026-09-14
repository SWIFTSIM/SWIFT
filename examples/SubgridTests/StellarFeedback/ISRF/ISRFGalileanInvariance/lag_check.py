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
gas (theory chapter, gas-frame paragraph of the Lagrangian form).

Measured lag, per band: the mass-weighted first moment of u over the whole
periodic box, relative to the star, along its direction of motion,
  lag = -sum_i m_i u_i xi_i / sum_i m_i u_i,  xi_i = minimum image of x_i - x_star,
with no radius window and no u > 0 filter. Prediction for a uniform c_hyp
and no dissipation: the steady moment balance of the discrete update,
  lag = v_rel tau a / (1 - exp(-a)),  tau = lambda / c_hyp,  a = dt / tau,
with lambda = 1/(kappa rho) from the code's dust opacity.

--gate (alpha_max = alpha_floor = 0): PASS if |lag/prediction - 1| <= 0.02 in
both bands at the last snapshot. --report-only (dissipation on): print the
displacement (lag - prediction)/h next to the 0.05 h reference, exit 0.
Any non-finite value: FAIL. The box must satisfy L >= 30 lambda in both bands
for the gate; a smaller box is reported as INVALID.
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
SIGMA_D_CGS = {"FUV": 9e-22, "LW": 1.5e-21}
MU_H = 1.4
M_H_CGS = 1.6726219e-24
# Grackle's SolarMetalFractionByMass default
GRACKLE_SOLAR_Z = 0.01295
REL_BAR = 0.02
H_REFERENCE = 0.05
MIN_BOX_OVER_LAMBDA = 30.0


def parse_options():
    parser = argparse.ArgumentParser(description=__doc__)
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
        star = f["/PartType4/Coordinates"][0, :].astype(np.float64)
    xi = pos - star
    xi -= L * np.round(xi / L)
    xi_par = xi @ v_hat
    out = dict(time=t, boxsize=L, h_med=float(np.median(h)), bands={})
    for band, sigma in SIGMA_D_CGS.items():
        kappa = sigma * (Z / GRACKLE_SOLAR_Z) / (MU_H * M_H_CGS)
        lam = 1.0 / (kappa * rho * um / ul**3) / ul
        finite = bool(np.all(np.isfinite(u[band])))
        M0 = float(np.sum(m * u[band]))
        lag = (
            float(-np.sum(m * u[band] * xi_par) / M0) if finite and M0 != 0 else np.nan
        )
        out["bands"][band] = dict(
            lag=lag,
            lambda_med=float(np.median(lam)),
            finite=finite,
            negative_mass_share=(
                float(np.sum(m * np.minimum(u[band], 0)) / M0) if M0 != 0 else np.nan
            ),
        )
    return out


def main():
    opt = parse_options()
    params = yaml.safe_load(open(os.path.join(opt.run, "used_parameters.yml")))
    fb = params["GEARFeedback"]
    c_hyp = float(fb["ISRF_c_hyp_pin_for_debugging"])
    if c_hyp <= 0.0:
        print("FAIL: this check needs ISRF_c_hyp_pin_for_debugging > 0 (uniform tau).")
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
    print(
        f"{opt.run}: v_rel={v_rel:.4g} km/s = {v_rel / c_hyp:.3g} c_hyp, c_hyp={c_hyp} km/s, "
        f"dt={dt:.4e}, alpha_max={fb['ISRF_dissipation_alpha_max']}, "
        f"alpha_floor={fb['ISRF_dissipation_alpha_floor']}"
    )

    verdict = "PASS"
    summary = dict(run=opt.run, v_rel_kms=v_rel, c_hyp_kms=c_hyp, dt=dt, snapshots={})
    for path in files[-3:]:
        res = lag_of(path, v_hat)
        last = path == files[-1]
        for band, r in res["bands"].items():
            tau = r["lambda_med"] / c_hyp
            a = dt / tau
            pred = v_rel * tau * a / (-np.expm1(-a))
            r.update(
                tau=tau,
                a=a,
                prediction=pred,
                t_over_tau=res["time"] / tau,
                box_over_lambda=res["boxsize"] / r["lambda_med"],
            )
            ratio = r["lag"] / pred
            disp_h = (r["lag"] - pred) / res["h_med"]
            r.update(ratio=ratio, displacement_h=disp_h)
            print(
                f"  t={res['time']:.3e} {band}: lag={r['lag']:.4e} pred={pred:.4e} "
                f"ratio-1={ratio - 1:+.4f} disp={disp_h:+.4f} h a={a:.4f} "
                f"t/tau={r['t_over_tau']:.1f} L/lambda={r['box_over_lambda']:.1f} "
                f"neg-share={r['negative_mass_share']:.3f}"
            )
            if not last:
                continue
            if not r["finite"] or not np.isfinite(ratio):
                verdict = "FAIL"
            elif opt.gate:
                if r["box_over_lambda"] < MIN_BOX_OVER_LAMBDA and verdict != "FAIL":
                    verdict = "INVALID"
                elif abs(ratio - 1.0) > REL_BAR:
                    verdict = "FAIL"
        summary["snapshots"][os.path.basename(path)] = res
    if opt.report_only and verdict == "PASS":
        verdict = "REPORT"
        print(f"  report only: displacement against the {H_REFERENCE} h reference")
    summary["verdict"] = verdict
    with open(os.path.join(opt.run, "lag_metrics.json"), "w") as f:
        json.dump(summary, f, indent=2, default=float)
    print(f"{opt.run}: lag check {verdict}")
    sys.exit(0 if verdict in ("PASS", "REPORT") else 1)


if __name__ == "__main__":
    main()
