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
Per-run metrics for the ISRFPropagationSpeedSweep example (c_hyp sweep, Legs
P/N/R/I -- see README). Computes the effective Courant number nu_eff from the
run's own measured h and dt (never assumed), the stability bound nu_max at
the run's own alpha, the pulse front position, and the realized-timestep
validity precondition Leg P's cross-run comparison depends on. Writes all of
it to --json-out; cross-run comparisons (self-similarity, retardation
collapse, the stability bracket) are done separately by sweep_compare.py.
"""

import argparse
import glob
import re
import sys

import h5py
import matplotlib
import numpy as np
import yaml

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SPEED_OF_LIGHT_KM_S = 2.99792458e5


def parse_options():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-s", "--snapshot", default="snap/snapshot_*.hdf5",
        help="Glob pattern for snapshots to consider (default: %(default)s)",
    )
    parser.add_argument(
        "--timesteps-log", default="timesteps.txt",
        help="Path to the run's own timesteps.txt (default: %(default)s).",
    )
    parser.add_argument(
        "--log", default="output.log",
        help="Path to the run's own stdout/stderr log (default: %(default)s).",
    )
    parser.add_argument(
        "--used-parameters", default="used_parameters.yml",
        help="Path to the run's own used_parameters.yml (default: %(default)s).",
    )
    parser.add_argument(
        "--c-hyp-pin", type=float, default=0.0,
        help="GEARFeedback:LW_FUV_c_hyp_pin_for_debugging used by the run, "
        "km/s (default: %(default)s; 0 = closure).",
    )
    parser.add_argument(
        "--c-hyp-margin", type=float, default=0.5,
        help="GEARFeedback:LW_FUV_c_hyp_margin used by the run (default: %(default)s).",
    )
    parser.add_argument("--n-bins", type=int, default=60, help="Radial bins.")
    parser.add_argument(
        "--expect-stable", dest="expect_stable", action="store_true", default=True,
        help="Gate M-P3: FAIL if the run destabilizes (default).",
    )
    parser.add_argument(
        "--expect-unstable", dest="expect_stable", action="store_false",
        help="Gate M-P3 the other way: destabilizing is CONFIRMED, not a failure "
        "(Leg N's N5/N6).",
    )
    parser.add_argument(
        "--output", default="isrf_propagation_speed_sweep_check.png",
        help="Output plot filename.",
    )
    parser.add_argument(
        "--json-out", default="sweep_metrics.json", help="Output metrics filename."
    )
    return parser.parse_args()


def modal_bulk_dt(path, n_total):
    """Modal interval between successive whole-box (Updates == n_total) steps;
    robust to a transient step-0 entry and to a run that never got that far."""
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
        return 0.0, 0
    full = times[updates == n_total]
    if len(full) < 2:
        vals, counts = np.unique(times[1:] - times[:-1], return_counts=True)
        dt = float(vals[np.argmax(counts)]) if len(vals) > 0 else 0.0
        return dt, int(times[-1] / dt) if dt > 0 else len(times)
    diffs = np.round(np.diff(full), 14)
    vals, counts = np.unique(diffs, return_counts=True)
    return float(vals[np.argmax(counts)]), len(full) - 1


def nu_max_of(alpha):
    """Positive root of 6.2*alpha*nu + 0.70*nu^2 = 2 (feedback_properties.h)."""
    a, b, c = 0.70, 6.2 * alpha, -2.0
    return (-b + np.sqrt(b * b - 4 * a * c)) / (2 * a)


def load_used_parameters(path):
    try:
        with open(path) as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        return {}


def outer_edge_above_threshold(r, u, threshold, n_bins, r_max):
    edges = np.linspace(0, r_max, n_bins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    idx = np.clip(np.digitize(r, edges) - 1, 0, n_bins - 1)
    means = np.full(n_bins, 0.0)
    for i in range(n_bins):
        sel = idx == i
        if sel.sum() > 0:
            means[i] = u[sel].mean()
    above = means > threshold
    return (float(centres[above].max()) if above.any() else 0.0), centres, means


def load_snapshot(path):
    with h5py.File(path, "r") as f:
        header = f["/Header"]
        time = float(np.asarray(header.attrs["Time"]).flat[0])
        boxsize = np.asarray(header.attrs["BoxSize"], dtype=float).flatten()[0]
        gas = f["/PartType0"]
        pos = gas["Coordinates"][:, :]
        h = gas["SmoothingLengths"][:].astype(np.float64)
        ids = gas["ParticleIDs"][:]
        rho = gas["Densities"][:].astype(np.float64)
        u_fuv = gas["FUVSpecificEnergies"][:].astype(np.float64)
        u_lw = gas["LWSpecificEnergies"][:].astype(np.float64)
    return dict(
        time=time, boxsize=boxsize, pos=pos, h=h, ids=ids, rho=rho,
        u_fuv=u_fuv, u_lw=u_lw,
    )


def main():
    opt = parse_options()
    files = sorted(glob.glob(opt.snapshot))
    if not files:
        raise RuntimeError(f"No snapshots match {opt.snapshot}")

    used_params = load_used_parameters(opt.used_parameters)
    fb = used_params.get("GEARFeedback", {}) if used_params else {}
    ti = used_params.get("TimeIntegration", {}) if used_params else {}
    alpha_max = float(fb.get("LW_FUV_dissipation_alpha_max", 0.25))
    alpha_pin = float(fb.get("LW_FUV_dissipation_alpha_pin_for_debugging", 0.0))
    dt_max_param = float(ti.get("dt_max", np.nan))
    alpha_eff = alpha_pin if alpha_pin > 0.0 else alpha_max

    snaps = [load_snapshot(fn) for fn in files]
    n_gas = snaps[0]["pos"].shape[0]

    dt_bulk, n_steps = modal_bulk_dt(opt.timesteps_log, n_gas)
    h_med_last = float(np.median(snaps[-1]["h"]))

    if opt.c_hyp_pin > 0.0:
        c_hyp = opt.c_hyp_pin
    else:
        c_hyp = min(opt.c_hyp_margin * h_med_last / dt_bulk, SPEED_OF_LIGHT_KM_S) \
            if dt_bulk > 0 else 0.0

    nu_eff = c_hyp * dt_bulk / h_med_last if h_med_last > 0 else 0.0
    nu_max = nu_max_of(alpha_eff)

    print(f"--- M-P1: nu_eff / stability bound ---")
    print(f"h_med (last snapshot) = {h_med_last:.6e}")
    print(f"dt_bulk (timesteps.txt, modal) = {dt_bulk:.6e}")
    print(f"c_hyp = {c_hyp:.6e} km/s")
    print(f"alpha_eff (pin if >0 else alpha_max) = {alpha_eff:.4f}")
    print(f"nu_eff = {nu_eff:.6f}")
    print(f"nu_max (bound at alpha_eff) = {nu_max:.6f}")
    print(f"nu_eff / nu_max = {nu_eff / nu_max if nu_max > 0 else float('nan'):.4f}")

    print(f"\n--- M-P6: realized-timestep validity (Leg P precondition) ---")
    print(f"dt_bulk = {dt_bulk:.6e}")
    print(f"realized step count = {n_steps}")
    print(f"dt_max (used_parameters.yml) = {dt_max_param:.6e}")
    print(f"dt_bulk / dt_max = {dt_bulk / dt_max_param if dt_max_param > 0 else float('nan'):.6f}")
    print(f"nu_eff = {nu_eff:.6f}")

    # M-P5: pin assertion.
    all_ok = True
    if opt.c_hyp_pin > 0.0:
        with open(opt.log) as f:
            log_text = f.read()
        if "LW_FUV_c_hyp_pin_for_debugging is set" not in log_text:
            all_ok = False
            print("\nFAIL M-P5: --c-hyp-pin > 0 but the pin warning is absent "
                  "from the log -- the pin did not take effect.")
        else:
            print("\nM-P5 PASS: pin warning present in the log.")

    # M-P2: gas-drift control. Snapshot particle order is not guaranteed
    # stable across outputs, so match by ParticleIDs before differencing.
    order0 = np.argsort(snaps[0]["ids"])
    orderN = np.argsort(snaps[-1]["ids"])
    t0_pos, t0_rho = snaps[0]["pos"][order0], snaps[0]["rho"][order0]
    tN_pos, tN_rho = snaps[-1]["pos"][orderN], snaps[-1]["rho"][orderN]
    boxsize = snaps[0]["boxsize"]
    dx = tN_pos - t0_pos
    dx -= boxsize * np.round(dx / boxsize)
    max_disp = float(np.sqrt(np.sum(dx**2, axis=1)).max())
    max_disp_h = max_disp / h_med_last if h_med_last > 0 else float("inf")
    max_drho = float(np.max(np.abs(tN_rho - t0_rho) / t0_rho))
    drift_void = max_disp_h > 0.1
    print(f"\n--- M-P2: gas-drift control ---")
    print(f"max displacement / h_med = {max_disp_h:.6f}"
          f"{'  -- VOID (>0.1h): invariance gates below do not apply' if drift_void else ''}")
    print(f"max |drho|/rho0 = {max_drho:.6f}")

    # M-P3: stability (max|u|, growth factor, NaNs).
    print(f"\n--- M-P3: stability ---")
    stability_ok = True
    for band, key in (("FUV", "u_fuv"), ("LW", "u_lw")):
        max_u_series = np.array([np.max(np.abs(s[key])) for s in snaps])
        n_nan = int(sum(np.sum(~np.isfinite(s[key])) for s in snaps))
        g = max_u_series[1:] / np.where(max_u_series[:-1] != 0, max_u_series[:-1], 1e-300)
        # 3 consecutive intervals (after the first) with g > 1.5.
        unstable = n_nan > 0
        if len(g) >= 4:
            for i in range(1, len(g) - 2):
                if g[i] > 1.5 and g[i + 1] > 1.5 and g[i + 2] > 1.5:
                    unstable = True
                    break
        print(f"{band}: max|u| per snapshot = {np.array2string(max_u_series, precision=3)}")
        print(f"{band}: growth factors = {np.array2string(g, precision=3)}, n_nan={n_nan}")
        if opt.expect_stable and unstable:
            stability_ok = False
            print(f"{band}: FAIL M-P3 -- run destabilized (NaN or sustained growth > 1.5x)")
        elif not opt.expect_stable and unstable:
            print(f"{band}: M-P3 CONFIRMED unstable, as expected for this run.")
        elif not opt.expect_stable and not unstable:
            print(f"{band}: run stayed stable; --expect-unstable was NOT confirmed.")
    all_ok &= stability_ok

    # M-P4: pulse front position.
    print(f"\n--- M-P4: pulse front position ---")
    last = snaps[-1]
    centre = np.array([0.5 * last["boxsize"]] * 3)
    dxc = last["pos"] - centre
    dxc -= last["boxsize"] * np.round(dxc / last["boxsize"])
    r = np.sqrt(np.sum(dxc**2, axis=1))
    r_max_plot = 0.45 * last["boxsize"]
    front_results = {}
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(snaps)))
    for band, key in (("FUV", "u_fuv"), ("LW", "u_lw")):
        u = last[key]
        u_max = float(u.max()) if u.size else 0.0
        band_results = {}
        for eps in (0.1, 0.01, 0.001):
            r_edge, centres, means = outer_edge_above_threshold(
                r, u, eps * u_max, opt.n_bins, r_max_plot
            )
            r_edge_h = r_edge / h_med_last if h_med_last > 0 else float("nan")
            denom = c_hyp * last["time"]
            r_edge_norm = r_edge / denom if denom > 0 else float("nan")
            band_results[eps] = dict(r_edge=r_edge, r_edge_h=r_edge_h, r_edge_norm=r_edge_norm)
            print(f"{band} eps={eps}: r_edge={r_edge:.4e} ({r_edge_h:.2f} h), "
                  f"r_edge/(c_hyp*t)={r_edge_norm:.4f}")
        front_results[band] = band_results

        ax = axes[0] if band == "FUV" else axes[1]
        for i, s in enumerate(snaps):
            dxi = s["pos"] - centre
            dxi -= s["boxsize"] * np.round(dxi / s["boxsize"])
            ri = np.sqrt(np.sum(dxi**2, axis=1))
            _, c_i, m_i = outer_edge_above_threshold(
                ri, s[key], 0.0, opt.n_bins, 0.45 * s["boxsize"]
            )
            valid = m_i > 0
            ax.semilogy(c_i[valid] / h_med_last, m_i[valid], "-",
                        color=colors[i], label=f"t={s['time']:.2e}")
        ax.axvline(c_hyp * last["time"] / h_med_last, color="k", ls="--", lw=1)
        ax.set_xlabel("r / h_med")
        ax.set_title(f"{band}: solid = u(r); dashed = c_hyp*t")
    axes[0].set_ylabel("u(r) (binned mean)")
    axes[0].legend(fontsize=6, loc="upper right")
    fig.tight_layout()
    fig.savefig(opt.output, dpi=150)
    print(f"\nPlot saved to {opt.output}")

    metrics = dict(
        h_med=h_med_last, dt_bulk=dt_bulk, n_steps=n_steps, dt_max=dt_max_param,
        c_hyp=c_hyp, c_hyp_pin=opt.c_hyp_pin, c_hyp_margin=opt.c_hyp_margin,
        alpha_max=alpha_max, alpha_pin=alpha_pin, alpha_eff=alpha_eff,
        nu_eff=nu_eff, nu_max=nu_max,
        max_disp_h=max_disp_h, max_drho=max_drho, drift_void=bool(drift_void),
        stability_ok=bool(stability_ok), expect_stable=bool(opt.expect_stable),
        front=front_results, time_last=last["time"],
    )
    import json
    with open(opt.json_out, "w") as f:
        json.dump(metrics, f, indent=2)
    print(f"Metrics written to {opt.json_out}")

    all_ok &= stability_ok
    if not all_ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
