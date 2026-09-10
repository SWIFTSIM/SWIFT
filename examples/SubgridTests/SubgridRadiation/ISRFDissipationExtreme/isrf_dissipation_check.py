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
Regression check for the LW/FUV Stage-1 triggered artificial-dissipation
term (`GEARFeedback:LW_FUV_dissipation_alpha_max`).

Adapted from the sibling `ISRFCausalReach` example's own check script (same
pinned-neighbour IC family and the same causal-reach metric); this version
adds the sign-closure metric that is this example's own point.

Two independent checks, both per snapshot and per band (`FUV`, `LW`):

1. Sign closure: `n_neg` (count of `FUVSpecificEnergies`/
   `LWSpecificEnergies` < 0, the artificially-heated particle included)
   and `u_min` are computed and reported every snapshot; the ratio
   `|u_min|/u_plateau` is gated at `--ratio-threshold` (2% by default) on
   snapshots in the run's last third only, matching the term's own design
   bar. This is the failure mode the dissipation term exists to remove.
2. Causal reach: the propagated field must not exceed its causal reach
   `r <= c_hyp*(t - t0)` (`t0` = the star's BirthTime, 0 here) by more than
   ordinary SPH kernel smearing. `c_hyp` is reconstructed from the run's own
   `timesteps.txt` and each snapshot's own median smoothing length, via the
   same formula `src/feedback/GEAR/radiation_isrf.c` uses:
   `c_hyp = min(LW_FUV_c_hyp_margin * h / dt, c)`. `C(eps)`, reported per
   epsilon in {0.1, 0.01, 0.001}, is how far past the causal front (in units
   of the local smoothing length) the field still exceeds
   `eps * u_plateau` -- an over-smoothing dissipation coefficient would
   inflate it.

`--hot-particle-id` (see `hot_particle_id.txt`, written by `makeIC.py`)
excludes the artificially-heated particle from the bulk h/dt estimate and
from `u_plateau` (its own field value is naturally far above the rest of
the box and would otherwise swamp both), and reports it separately.
"""

import argparse
import glob
import sys

import h5py
import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

SPEED_OF_LIGHT_KM_S = 2.99792458e5


def parse_options():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-s",
        "--snapshot",
        default="snap/snapshot_*.hdf5",
        help="Glob pattern for snapshots to consider (default: %(default)s)",
    )
    parser.add_argument(
        "--timesteps-log",
        default="timesteps.txt",
        help="Path to the run's own timesteps.txt (default: %(default)s).",
    )
    parser.add_argument(
        "--c-hyp-margin",
        type=float,
        default=0.5,
        help="GEARFeedback:LW_FUV_c_hyp_margin used by the run (default: %(default)s).",
    )
    parser.add_argument(
        "--hot-particle-id",
        type=int,
        default=-1,
        help="Pinned-neighbour variant: ID of the artificially-heated gas "
        "particle (see hot_particle_id.txt), excluded from the bulk h/dt "
        "and u_plateau estimates and reported on separately.",
    )
    parser.add_argument("--n-bins", type=int, default=120, help="Radial bins.")
    parser.add_argument(
        "--near-source-h",
        type=float,
        default=4.5,
        help="Exclusion radius (in units of h) around the star for the "
        "monotonicity check, to skip ordinary injection-kernel ripples "
        "(default: %(default)s).",
    )
    parser.add_argument(
        "--fail-margin-h",
        type=float,
        default=20.0,
        help="Hard-fail threshold (in units of h) beyond c_hyp*t: any level "
        "of u beyond this is a real finding, not smearing.",
    )
    parser.add_argument(
        "--ratio-threshold",
        type=float,
        default=0.02,
        help="Sign-closure pass bar: |u_min|/u_plateau, gated only on "
        "snapshots in the run's last third (default: %(default)s, i.e. 2%%).",
    )
    parser.add_argument(
        "--output",
        default="isrf_dissipation_check.png",
        help="Output plot filename.",
    )
    return parser.parse_args()


def modal_bulk_dt(path, n_total):
    """The bulk gas's own physical time-step, as the modal interval between
    successive whole-box updates (Updates == n_total): robust even when a
    single pinned particle's own much finer ticks dominate the plain modal
    Time-step column (the pinned-neighbour variant's whole point)."""
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
    full = times[updates == n_total]
    if len(full) < 2:
        # No pinned particle: every active step already updates everyone.
        vals, counts = np.unique(times[1:] - times[:-1], return_counts=True)
        return float(vals[np.argmax(counts)])
    diffs = np.round(np.diff(full), 12)
    vals, counts = np.unique(diffs, return_counts=True)
    return float(vals[np.argmax(counts)])


def load_snapshot(path):
    with h5py.File(path, "r") as f:
        header = f["/Header"]
        time = float(np.asarray(header.attrs["Time"]).flat[0])
        boxsize = np.asarray(header.attrs["BoxSize"], dtype=float).flatten()[0]
        gas = f["/PartType0"]
        pos = gas["Coordinates"][:, :]
        h = gas["SmoothingLengths"][:].astype(np.float64)
        ids = gas["ParticleIDs"][:]
        u_fuv = gas["FUVSpecificEnergies"][:].astype(np.float64)
        u_lw = gas["LWSpecificEnergies"][:].astype(np.float64)
        star = f["/PartType4"]
        star_pos = star["Coordinates"][0, :]
    return dict(
        time=time,
        boxsize=boxsize,
        pos=pos,
        h=h,
        ids=ids,
        u_fuv=u_fuv,
        u_lw=u_lw,
        star_pos=star_pos,
    )


def radial_distance(pos, star_pos, boxsize):
    dx = pos - star_pos
    dx -= boxsize * np.round(dx / boxsize)
    return np.sqrt(np.sum(dx**2, axis=1))


def outer_edge_above_threshold(r, u, threshold, n_bins, r_max):
    """Largest bin-centre radius where the radially-binned mean u still
    exceeds `threshold`; 0 if no bin exceeds it."""
    edges = np.linspace(0, r_max, n_bins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    idx = np.digitize(r, edges) - 1
    idx = np.clip(idx, 0, n_bins - 1)
    means = np.full(n_bins, 0.0)
    for i in range(n_bins):
        sel = idx == i
        if sel.sum() > 0:
            means[i] = u[sel].mean()
    above = means > threshold
    return float(centres[above].max()) if above.any() else 0.0, centres, means


def check_band(
    band, r, u, u_all_incl_hot, h_med, c_hyp, t, t0, opt, r_max_plot, in_last_third
):
    """`u`/`r` are the masked (bulk) arrays used for u_plateau and the
    causal-reach metric; `u_all_incl_hot` is the full, unmasked array used
    for sign closure -- a negative value at the heated particle itself is
    still a failure. `n_neg` is always reported; only the ratio is gated,
    and only on the run's last third (the pass bar this dissipation term is
    designed against: ratio <= --ratio-threshold there, in both bands)."""
    u_plateau = float(u.max())
    ok = True
    findings = []

    n_neg = int(np.sum(u_all_incl_hot < 0))
    u_min = float(u_all_incl_hot.min())
    ratio = abs(u_min) / u_plateau if u_plateau > 0 and u_min < 0 else 0.0
    if in_last_third and ratio > opt.ratio_threshold:
        ok = False
        findings.append(
            f"{band}: FAIL sign closure -- ratio |u_min|/u_plateau = "
            f"{ratio*100:.3f}% exceeds {opt.ratio_threshold*100:.1f}% "
            f"({n_neg} particles negative, u_min = {u_min:.3e})"
        )

    r_front = c_hyp * max(t - t0, 0.0)
    results = {}
    for eps in (0.1, 0.01, 0.001):
        r_edge, centres, means = outer_edge_above_threshold(
            r, u, eps * u_plateau, opt.n_bins, r_max_plot
        )
        C_h = (r_edge - r_front) / h_med if h_med > 0 else np.inf
        results[eps] = (r_edge, C_h)
        if r_edge > r_front + opt.fail_margin_h * h_med:
            ok = False
            findings.append(
                f"{band}: FAIL causal reach at epsilon={eps} -- edge at "
                f"{r_edge:.4e} ({(r_edge)/h_med:.2f} h) exceeds "
                f"c_hyp*t + {opt.fail_margin_h}h "
                f"({r_front + opt.fail_margin_h * h_med:.4e})"
            )

    # Monotonicity behind the front, outside the near-source exclusion zone.
    edges = np.linspace(0, r_max_plot, opt.n_bins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    idx = np.clip(np.digitize(r, edges) - 1, 0, opt.n_bins - 1)
    means = np.full(opt.n_bins, np.nan)
    for i in range(opt.n_bins):
        sel = idx == i
        if sel.sum() > 0:
            means[i] = u[sel].mean()
    exclude = centres < opt.near_source_h * h_med
    keep = ~np.isnan(means) & ~exclude
    worst_bump = 0.0
    worst_bump_r = None
    if keep.sum() > 1:
        vals = means[keep]
        rs = centres[keep]
        rel_increase = np.diff(vals) / np.maximum(np.abs(vals[:-1]), 1e-300)
        if len(rel_increase) > 0:
            j = int(np.argmax(rel_increase))
            worst_bump = float(rel_increase[j])
            worst_bump_r = float(rs[j])

    return dict(
        ok=ok,
        findings=findings,
        u_plateau=u_plateau,
        r_front=r_front,
        n_neg=n_neg,
        u_min=u_min,
        ratio=ratio,
        results=results,
        worst_bump=worst_bump,
        worst_bump_r=worst_bump_r,
        centres=centres,
        means=means,
    )


def main():
    opt = parse_options()
    files = sorted(glob.glob(opt.snapshot))
    if not files:
        raise RuntimeError(f"No snapshots match {opt.snapshot}")

    with h5py.File(files[0], "r") as f:
        n_gas = f["/PartType0/Coordinates"].shape[0]
    # "Updates" in timesteps.txt counts gas particles only (not the star).
    dt_bulk = modal_bulk_dt(opt.timesteps_log, n_gas)
    print(f"Bulk gas time-step (timesteps.txt, n_gas={n_gas}): {dt_bulk:.6e}")

    all_ok = True
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(files)))

    t0 = 0.0  # star BirthTime = 0 in this example's ICs (see README).

    # Last-third window for the gated sign-closure ratio: found from the
    # snapshots' own times, not assumed from TimeIntegration:time_end.
    times_all = [load_snapshot(fn)["time"] for fn in files]
    t_max = max(times_all)
    last_third_start = t0 + (2.0 / 3.0) * (t_max - t0)
    print(
        f"Last-third window (gated sign-closure ratio): "
        f"t >= {last_third_start:.4e} (t_max = {t_max:.4e})"
    )

    for i, fn in enumerate(files):
        snap = load_snapshot(fn)
        t = snap["time"]
        if t <= t0:
            continue
        in_last_third = t >= last_third_start
        pos, h, ids = snap["pos"], snap["h"], snap["ids"]
        gas_mask = np.ones(len(ids), dtype=bool)
        if opt.hot_particle_id >= 0:
            gas_mask = ids != opt.hot_particle_id
        h_med = float(np.median(h[gas_mask]))
        c_hyp = min(opt.c_hyp_margin * h_med / dt_bulk, SPEED_OF_LIGHT_KM_S)

        r_all = radial_distance(pos, snap["star_pos"], snap["boxsize"])
        r_max_plot = 0.45 * snap["boxsize"]

        print(
            f"\n--- t={t:.4e}  h_med={h_med:.4e}  c_hyp={c_hyp:.4f}  "
            f"r_front={c_hyp * (t - t0):.4e} ({c_hyp * (t - t0) / h_med:.2f} h)  "
            f"last_third={in_last_third} ---"
        )

        for band, u_field in (("FUV", "u_fuv"), ("LW", "u_lw")):
            u_all = snap[u_field]
            r, u = r_all[gas_mask], u_all[gas_mask]
            res = check_band(
                band, r, u, u_all, h_med, c_hyp, t, t0, opt, r_max_plot, in_last_third
            )
            all_ok &= res["ok"]
            status = "PASS" if res["ok"] else "FAIL"
            eps_str = ", ".join(
                f"eps={e}: edge={res['results'][e][0]/h_med:.2f}h "
                f"(C={res['results'][e][1]:.2f})"
                for e in (0.1, 0.01, 0.001)
            )
            print(
                f"{band}: n_neg={res['n_neg']}  u_min={res['u_min']:.4e}  "
                f"ratio={res['ratio']*100:.3f}%  u_plateau={res['u_plateau']:.4e}  "
                f"{eps_str}  worst_bump={res['worst_bump']:.3f} at "
                f"r={res['worst_bump_r']}  -> {status}"
            )
            for finding in res["findings"]:
                print("  " + finding)

            ax = axes[0] if band == "FUV" else axes[1]
            valid = res["means"] > 0
            ax.semilogy(
                res["centres"][valid] / h_med,
                res["means"][valid],
                "-",
                color=colors[i],
                label=f"t={t:.2e}",
            )
            ax.axvline(res["r_front"] / h_med, color=colors[i], ls="--", lw=1)

    for ax, band in zip(axes, ("FUV", "LW")):
        ax.set_xlabel("r / h")
        ax.set_title(f"{band}: solid = u(r); dashed = c_hyp*(t-t0)")
    axes[0].set_ylabel("u(r) (binned mean)")
    axes[0].legend(fontsize=6, loc="upper right")
    fig.tight_layout()
    fig.savefig(opt.output, dpi=150)
    print(f"\nPlot saved to {opt.output}")

    if not all_ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
