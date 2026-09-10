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
Per-run metrics for ISRFShearAsymmetry: the reflection-symmetry asymmetry
(A_centroid, A_spread, A_energy), the x-structure metric (slab geometry),
and the KH-contamination validity control. Pass/fail against the zero-shear
control is done separately by shear_compare.py, since the gate is
self-calibrating (2.8) and needs the sibling v_rel=0 run's own metrics.

A_centroid and A_spread are evaluated on the positive part of the field,
`w = m*max(u, 0)`, not the raw `w = m*u`: a large negative-weight share
makes the raw second moment an ill-conditioned difference of two similar
numbers (observed directly: a 31%-negative-weight run collapsed A_spread
to an outlier). The raw weighting is still computed and reported
alongside (`*_raw` fields, and printed) for comparison. A run whose
negative weight share exceeds NEG_WEIGHT_VOID_THRESHOLD in either blob
is marked void, folded into the same top-level `void` flag the
KH-contamination check uses, so `shear_compare.py`'s existing gate skips
it without further changes. A_energy keeps the raw weight: it is what
conservation acts on.
"""

import argparse
import glob
import json
import sys

import h5py
import numpy as np
import yaml

SPEED_OF_LIGHT_KM_S = 2.99792458e5
NEG_WEIGHT_VOID_THRESHOLD = 0.10  # negative-weight share above which A_centroid/A_spread are voided


def parse_options():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-s", "--snapshot", default="snap/snapshot_*.hdf5")
    parser.add_argument("--timesteps-log", default="timesteps.txt")
    parser.add_argument("--used-parameters", default="used_parameters.yml")
    parser.add_argument(
        "--variant", choices=["shear", "contrast"], default="shear",
    )
    parser.add_argument(
        "--source-geometry", choices=["blobs", "slab"], default="blobs",
    )
    parser.add_argument("--v-shear", type=float, default=1.0, help="km/s")
    parser.add_argument("--c-hyp-pin", type=float, default=0.0)
    parser.add_argument("--c-hyp-margin", type=float, default=0.5)
    parser.add_argument("--pulse-sigma-h", type=float, default=2.0)
    parser.add_argument("--n-x-bins", type=int, default=16)
    parser.add_argument("--n-y-slabs", type=int, default=32)
    parser.add_argument("--json-out", default="shear_metrics.json")
    return parser.parse_args()


def modal_bulk_dt(path, n_total):
    times, updates = [], []
    with open(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split()
            step = int(parts[0])
            if step == 0:
                continue
            times.append(float(parts[1]))
            updates.append(int(parts[7]))
    times, updates = np.array(times), np.array(updates)
    if len(times) == 0:
        return 0.0
    full = times[updates == n_total]
    if len(full) < 2:
        vals, counts = np.unique(times[1:] - times[:-1], return_counts=True)
        return float(vals[np.argmax(counts)]) if len(vals) else 0.0
    diffs = np.round(np.diff(full), 14)
    vals, counts = np.unique(diffs, return_counts=True)
    return float(vals[np.argmax(counts)])


def load_snapshot(path):
    with h5py.File(path, "r") as f:
        header = f["/Header"]
        time = float(np.asarray(header.attrs["Time"]).flat[0])
        boxsize = np.asarray(header.attrs["BoxSize"], dtype=float).flatten()[0]
        gas = f["/PartType0"]
        pos = gas["Coordinates"][:, :]
        vel = gas["Velocities"][:, :]
        mass = gas["Masses"][:].astype(np.float64)
        rho = gas["Densities"][:].astype(np.float64)
        h = gas["SmoothingLengths"][:].astype(np.float64)
        ids = gas["ParticleIDs"][:]
        u_int = gas["InternalEnergies"][:].astype(np.float64)
        u_fuv = gas["FUVSpecificEnergies"][:].astype(np.float64)
        u_lw = gas["LWSpecificEnergies"][:].astype(np.float64)
    return dict(
        time=time, boxsize=boxsize, pos=pos, vel=vel, mass=mass, rho=rho, h=h,
        ids=ids, u_int=u_int, u_fuv=u_fuv, u_lw=u_lw,
    )


def min_image(dx, boxsize):
    return dx - boxsize * np.round(dx / boxsize)


def blob_membership(pos0, boxsize, L_code, h_mean, sigma_h):
    """Initial-position membership (fixed for the whole run), matching
    makeIC.py's blob centres and 3*sigma radius exactly."""
    sigma = sigma_h * h_mean
    centre_A = np.array([0.25 * L_code, 0.50 * L_code, 0.50 * L_code])
    centre_B = np.array([0.75 * L_code, 0.00 * L_code, 0.50 * L_code])
    dA = min_image(pos0 - centre_A, boxsize)
    dB = min_image(pos0 - centre_B, boxsize)
    rA = np.sqrt(np.sum(dA**2, axis=1))
    rB = np.sqrt(np.sum(dB**2, axis=1))
    return (rA <= 3.0 * sigma), (rB <= 3.0 * sigma), sigma, centre_A, centre_B


def main():
    opt = parse_options()
    files = sorted(glob.glob(opt.snapshot))
    if not files:
        raise RuntimeError(f"No snapshots match {opt.snapshot}")

    with open(opt.used_parameters) as f:
        used = yaml.safe_load(f)
    fb = used.get("GEARFeedback", {})
    ti = used.get("TimeIntegration", {})
    dt_max_param = float(ti.get("dt_max", np.nan))
    alpha_max = float(fb.get("LW_FUV_dissipation_alpha_max", 0.25))
    alpha_pin = float(fb.get("LW_FUV_dissipation_alpha_pin_for_debugging", 0.0))

    snap0 = load_snapshot(files[0])
    snap_last = load_snapshot(files[-1])
    n_gas = snap0["pos"].shape[0]
    L_code = snap0["boxsize"]
    h_mean_analytic = 1.2348 * L_code / n_gas ** (1.0 / 3.0)
    h_med_last = float(np.median(snap_last["h"]))

    dt_bulk = modal_bulk_dt(opt.timesteps_log, n_gas)
    if opt.c_hyp_pin > 0.0:
        c_hyp = opt.c_hyp_pin
    else:
        c_hyp = min(opt.c_hyp_margin * h_med_last / dt_bulk, SPEED_OF_LIGHT_KM_S) \
            if dt_bulk > 0 else 0.0
    nu_eff = c_hyp * dt_bulk / h_med_last if h_med_last > 0 else 0.0

    # M-S7 configuration echo.
    c_s = float(np.sqrt(5.0 / 3.0 * 2.0 / 3.0 * snap0["u_int"].mean()))  # gamma=5/3, u=kB*T/((g-1)*mu*mp) -> c_s^2=(g)(g-1)u
    mach = opt.v_shear / c_s if c_s > 0 else float("nan")
    print("--- M-S7: configuration echo ---")
    print(f"h_med (last) = {h_med_last:.6e}   h_mean (analytic, makeIC.py's own) = {h_mean_analytic:.6e}")
    print(f"dt_bulk = {dt_bulk:.6e}   dt_max (param) = {dt_max_param:.6e}")
    print(f"c_hyp = {c_hyp:.6e} km/s   nu_eff = {nu_eff:.6f}")
    print(f"c_s (from mean InternalEnergy) = {c_s:.4f} km/s   Mach = v_shear/c_s = {mach:.4f}")

    order0 = np.argsort(snap0["ids"])
    orderL = np.argsort(snap_last["ids"])
    ids0 = snap0["ids"][order0]
    pos0 = snap0["pos"][order0]
    mass0 = snap0["mass"][order0]
    rho0 = snap0["rho"][order0]

    idsL = snap_last["ids"][orderL]
    if not np.array_equal(ids0, idsL):
        raise RuntimeError("ParticleIDs mismatch between first and last snapshot.")

    # M-S6: KH-contamination control, per snapshot (report the last one).
    print("\n--- M-S6: KH-contamination control ---")
    kh_last, drho_last = None, None
    for fn in files:
        s = load_snapshot(fn)
        order = np.argsort(s["ids"])
        rho_s = s["rho"][order]
        vy_rms = float(np.sqrt(np.mean(s["vel"][:, 1] ** 2)))
        if opt.v_shear > 0:
            kh = vy_rms / max(abs(opt.v_shear), 1e-30)
            kh_str = f"kh={kh:.4f}"
        else:
            kh = vy_rms / c_s if c_s > 0 else float("nan")
            kh_str = f"rms(v_y)/c_s={kh:.4f} (v_shear=0, kh undefined)"
        drho = float(np.max(np.abs(rho_s - rho0) / rho0))
        print(f"t={s['time']:.4e}: {kh_str}  drho={drho:.4f}")
        kh_last, drho_last = kh, drho
    void = (opt.v_shear > 0 and kh_last is not None and kh_last > 0.05) or (
        drho_last is not None and drho_last > 0.10
    )
    print(f"Validity: {'VOID' if void else 'OK'} (kh>0.05 or drho>0.10 at last snapshot)")

    # M-S0: population imbalance (report; gate applied once at the smoke-test stage).
    metrics = dict(
        h_med=h_med_last, h_mean_analytic=h_mean_analytic, dt_bulk=dt_bulk,
        dt_max=dt_max_param, c_hyp=c_hyp, nu_eff=nu_eff, c_s=c_s, mach=mach,
        alpha_max=alpha_max, alpha_pin=alpha_pin,
        kh_last=kh_last, drho_last=drho_last, void=bool(void),
        time_last=snap_last["time"], v_shear=opt.v_shear,
        # "positive-part": A_centroid/A_spread use w=m*max(u,0); a JSON file
        # without this field predates the fix and used raw w=m*u instead --
        # regenerate it before gating against a run using this version.
        moment_weighting="positive-part",
    )

    if opt.source_geometry == "blobs":
        in_A, in_B, sigma, centre_A, centre_B = blob_membership(
            pos0, L_code, L_code, h_mean_analytic, opt.pulse_sigma_h
        )
        N_A, N_B = int(in_A.sum()), int(in_B.sum())
        print(f"\n--- M-S0: population imbalance (report; gated once at smoke-test) ---")
        print(f"sigma (script) = {sigma:.6e}, centre_A={centre_A}, centre_B={centre_B}")
        print(f"N_A={N_A}  N_B={N_B}  (N_A-N_B)/(N_A+N_B)={(N_A - N_B) / (N_A + N_B):.4%}")

        # Comoving-frame centroid/spread/energy asymmetry (M-S1/S2/S3).
        s_A, s_B = +1.0, -1.0
        u_last_fuv = snap_last["u_fuv"][orderL]
        u_last_lw = snap_last["u_lw"][orderL]
        pos_last = snap_last["pos"][orderL]
        t = snap_last["time"]

        def moments(w, dx):
            w_sum = w.sum()
            if w_sum == 0:
                return 0.0, 0.0
            Dx = float(np.sum(w * dx[:, 0]) / w_sum)
            Sxx = float(np.sum(w * dx[:, 0] ** 2) / w_sum)
            return Dx, Sxx

        results = {}
        void_neg_weight_any = False
        for band, u_last in (("FUV", u_last_fuv), ("LW", u_last_lw)):
            band_res = {}
            for label, in_blob, s_b, centre in (
                ("A", in_A, s_A, centre_A), ("B", in_B, s_B, centre_B)
            ):
                x_ref = centre + np.array([s_b * (opt.v_shear / 2.0) * t, 0.0, 0.0])
                dx = min_image(pos_last[in_blob] - x_ref, L_code)
                w_raw = mass0[in_blob] * u_last[in_blob]
                w_clip = mass0[in_blob] * np.maximum(u_last[in_blob], 0.0)
                w_abs_sum = float(np.sum(np.abs(w_raw)))
                neg_share = float(-np.sum(w_raw[w_raw < 0]) / w_abs_sum) if w_abs_sum > 0 else 0.0
                Dx, Sxx = moments(w_clip, dx)
                Dx_raw, Sxx_raw = moments(w_raw, dx)
                E = float(w_raw.sum())
                band_res[label] = dict(
                    Dx=Dx, Sxx=Sxx, E=E,
                    Dx_raw=Dx_raw, Sxx_raw=Sxx_raw,
                    neg_weight_share=neg_share,
                )
            Dx_A, Dx_B = band_res["A"]["Dx"], band_res["B"]["Dx"]
            Sxx_A, Sxx_B = band_res["A"]["Sxx"], band_res["B"]["Sxx"]
            E_A, E_B = band_res["A"]["E"], band_res["B"]["E"]
            Dx_A_raw, Dx_B_raw = band_res["A"]["Dx_raw"], band_res["B"]["Dx_raw"]
            Sxx_A_raw, Sxx_B_raw = band_res["A"]["Sxx_raw"], band_res["B"]["Sxx_raw"]
            neg_A = band_res["A"]["neg_weight_share"]
            neg_B = band_res["B"]["neg_weight_share"]
            void_neg_weight = max(neg_A, neg_B) > NEG_WEIGHT_VOID_THRESHOLD
            void_neg_weight_any |= void_neg_weight

            A_centroid = (Dx_A + Dx_B) / h_med_last if h_med_last > 0 else float("nan")
            A_spread = 2 * (Sxx_A - Sxx_B) / (Sxx_A + Sxx_B) if (Sxx_A + Sxx_B) else float("nan")
            A_energy = 2 * (E_A - E_B) / (E_A + E_B) if (E_A + E_B) else float("nan")
            A_centroid_raw = (Dx_A_raw + Dx_B_raw) / h_med_last if h_med_last > 0 else float("nan")
            A_spread_raw = (
                2 * (Sxx_A_raw - Sxx_B_raw) / (Sxx_A_raw + Sxx_B_raw)
                if (Sxx_A_raw + Sxx_B_raw) else float("nan")
            )
            print(f"{band}: neg-weight share A={neg_A:.1%} B={neg_B:.1%}"
                  f"  ({'VOID, >' if void_neg_weight else 'OK, <='}"
                  f"{NEG_WEIGHT_VOID_THRESHOLD:.0%})")
            print(f"{band}: Dx_A={Dx_A:.4e} Dx_B={Dx_B:.4e} -> A_centroid={A_centroid:.4e}"
                  f"  (raw weighting: {A_centroid_raw:.4e})")
            print(f"{band}: Sxx_A={Sxx_A:.4e} Sxx_B={Sxx_B:.4e} -> A_spread={A_spread:.4e}"
                  f"  (raw weighting: {A_spread_raw:.4e})")
            print(f"{band}: E_A={E_A:.4e} E_B={E_B:.4e} -> A_energy={A_energy:.4e}")
            results[band] = dict(
                Dx_A=Dx_A, Dx_B=Dx_B, A_centroid=A_centroid,
                Sxx_A=Sxx_A, Sxx_B=Sxx_B, A_spread=A_spread,
                E_A=E_A, E_B=E_B, A_energy=A_energy,
                A_centroid_raw=A_centroid_raw, A_spread_raw=A_spread_raw,
                neg_weight_share_A=neg_A, neg_weight_share_B=neg_B,
                void_neg_weight=bool(void_neg_weight),
            )
        metrics["N_A"] = N_A
        metrics["N_B"] = N_B
        metrics["bands"] = results
        metrics["void_neg_weight"] = bool(void_neg_weight_any)
        void = void or void_neg_weight_any
        metrics["void"] = bool(void)

        # Total conservation control (report only).
        for band, u_field in (("FUV", "u_fuv"), ("LW", "u_lw")):
            tot0 = float(np.sum(mass0 * snap0[u_field][order0]))
            totL = float(np.sum(mass0 * snap_last[u_field][orderL]))
            drift = (totL - tot0) / tot0 if tot0 != 0 else float("nan")
            print(f"{band}: total sum(m*u) t0={tot0:.6e} tN={totL:.6e} drift={drift:.4%}")
            metrics.setdefault("total_drift", {})[band] = drift

    else:
        # M-S4: A_xstructure (slab geometry, report only).
        pos_last = snap_last["pos"][orderL]
        results = {}
        for band, key in (("FUV", "u_fuv"), ("LW", "u_lw")):
            u_last = snap_last[key][orderL]
            y_edges = np.linspace(0, L_code, opt.n_y_slabs + 1)
            x_edges = np.linspace(0, L_code, opt.n_x_bins + 1)
            y_idx = np.clip(np.digitize(pos_last[:, 1], y_edges) - 1, 0, opt.n_y_slabs - 1)
            worst = 0.0
            for iy in range(opt.n_y_slabs):
                sel = y_idx == iy
                if sel.sum() < opt.n_x_bins:
                    continue
                x_idx = np.clip(np.digitize(pos_last[sel, 0], x_edges) - 1, 0, opt.n_x_bins - 1)
                means = np.full(opt.n_x_bins, np.nan)
                for ix in range(opt.n_x_bins):
                    m = x_idx == ix
                    if m.sum() > 0:
                        means[ix] = u_last[sel][m].mean()
                valid = ~np.isnan(means)
                if valid.sum() > 1 and np.mean(means[valid]) != 0:
                    ratio = np.std(means[valid]) / abs(np.mean(means[valid]))
                    worst = max(worst, ratio)
            print(f"{band}: A_xstructure = {worst:.4e}")
            results[band] = dict(A_xstructure=worst)
        metrics["bands"] = results

    with open(opt.json_out, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"\nMetrics written to {opt.json_out}")

    if void:
        reasons = []
        if kh_last is not None and (
            (opt.v_shear > 0 and kh_last > 0.05) or (drho_last is not None and drho_last > 0.10)
        ):
            reasons.append("KH contamination")
        if opt.source_geometry == "blobs" and void_neg_weight_any:
            reasons.append(f"negative-weight share > {NEG_WEIGHT_VOID_THRESHOLD:.0%}")
        print(f"\nVALIDITY WARNING: this run is VOID ({', '.join(reasons)}) -- "
              "M-S1/S2/S3 results above should not be trusted.")


if __name__ == "__main__":
    main()
