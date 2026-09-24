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
Sec 6.2's causal-reach test: the propagated PE/LW field must not exceed its
causal reach `r <= c_hyp*t` (from injection start) beyond the SPH kernel's
own expected numerical smearing.

`c_hyp` is not itself a snapshot field (this project's "no ad hoc debug
macros" convention rules out adding one purely for this check).
It is instead reconstructed from the run's own `timesteps.txt` (the modal
gas time-step, excluding the transient step 0) and each snapshot's own
median smoothing length, via the exact formula radiation_isrf.c uses:
`c_hyp = min(ISRF_c_hyp_margin * h / dt, c)`. In the standard (uniform)
variant every gas particle shares (to float precision) the same h and dt, so
one scalar `c_hyp` applies to the whole box; the pinned-neighbour variant
uses `--hot-particle-id` to exclude that one particle from the "bulk" h/dt
estimate and from every metric below (see `--hot-particle-id`'s own help).

Pass criterion (Sec 6.2): for each of epsilon in {0.1, 0.01, 0.001}, the
outer radius where the radially-binned u(r) last exceeds
`epsilon * u_plateau` (u_plateau := max(u) this snapshot, no 1-D-style
flat plateau exists for a 3-D point source, so this is the most permissive
definition, per this project's own review of the 1-D reference numbers)
must not exceed `c_hyp*(t-t0) + C(epsilon)*h` by more than a generous
margin; C(epsilon) is measured and reported, not assumed, since the 1-D
reference values (this project's 1-D timestepping-stability verification,
Part F.2, C~1/6/12 h at 10%/1%/0.1%) do not necessarily transfer to a 3-D
point source's geometric (1/r) dilution. `t0` is the star's BirthTime (0
here, so no offset is applied): the run's own README/params.yml choice is
recorded, not assumed universally correct.

Negativity is GATED for this leg (unlike the Tier-1 steady-state check's
"informational only"): Sec 6.2 explicitly names the wavefront region as
where Sec 2.3's estimator weakness is most likely to actually show up.

`worst_bump`, the largest bin-to-bin relative rise of u(r) outside a
near-injection-kernel exclusion zone, is REPORTED ONLY and gates nothing.
It is a ratio of consecutive tail-bin means, and those means fall towards
zero away from the source, so its own noise is unbounded and no mechanism
supplies a bar for it. A bar taken from its observed spread would be a bar
fitted to the measurement.

Radial binning EXCLUDES particles outside the binned range rather than
folding them into the outermost bin. The range is the full periodic
minimum-image reach, `sqrt(3)/2 * L`, so no part of the box is outside it
and the causal-reach gate is never blind to a far-field violation. Folding
the box corners into the outermost bin would drag that bin's mean down with
their near-zero u and make the gate harder to trip the further out the
violation sits, which is the opposite of what this check is for.

Every field a gate reads (each band's specific energy, smoothing lengths,
positions, box size, time, the bulk time-step) is checked for finiteness
before it reaches a gate. A NaN compares False against every threshold, so
an unchecked NaN would silently turn both the negativity gate and the
causal-reach gate into a trivial pass; this check fails loudly instead,
naming the snapshot, the band and the count of bad values.
"""

import argparse
import glob
import re
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
        help="GEARFeedback:ISRF_c_hyp_margin used by the run (default: %(default)s).",
    )
    parser.add_argument(
        "--hot-particle-id",
        type=int,
        default=-1,
        help="Pinned-neighbour variant: ID of the artificially-heated gas "
        "particle (see hot_particle_id.txt), excluded from the bulk h/dt "
        "estimate and from every metric, including the negativity gate.",
    )
    parser.add_argument(
        "--n-bins",
        type=int,
        default=231,
        help="Radial bins over the full periodic reach sqrt(3)/2 L. The "
        "default keeps the same bin WIDTH as 120 bins over 0.45 L did, so "
        "the front position is quantised as before (default: %(default)s).",
    )
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
        "of u beyond this is a real finding, not smearing (Sec 6.2).",
    )
    parser.add_argument(
        "--output",
        default="isrf_causal_reach_check.png",
        help="Output plot filename.",
    )
    return parser.parse_args()


def check_finite(name, arr, snapshot, band=None):
    """Fail loudly if `arr` contains a NaN or an infinite value.

    A pass/fail check must never let a non-finite value slip through
    unnoticed: NaN compares False against every threshold, so an
    undetected NaN silently turns every downstream gate that reads it
    into a trivial pass rather than a real failure.

    Parameters
    ----------
    name : str
        Human-readable name of the quantity being checked.
    arr : numpy.ndarray
        Values to test for finiteness.
    snapshot : str
        Path of the snapshot the values came from, for the error message.
    band : str, optional
        Radiation band the quantity belongs to, if any.

    Raises
    ------
    RuntimeError
        If any element of `arr` is not finite.
    """
    arr = np.asarray(arr)
    n_bad = int(np.sum(~np.isfinite(arr)))
    if n_bad > 0:
        where = f", band {band}" if band else ""
        raise RuntimeError(
            f"{snapshot}{where}: {n_bad}/{arr.size} non-finite value(s) "
            f"(NaN or inf) in {name} -- refusing to gate on a corrupted "
            f"field."
        )


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
        u_pe = gas["PESpecificEnergies"][:].astype(np.float64)
        u_lw = gas["LWSpecificEnergies"][:].astype(np.float64)
        star = f["/PartType4"]
        star_pos = star["Coordinates"][0, :]
    return dict(
        time=time,
        boxsize=boxsize,
        pos=pos,
        h=h,
        ids=ids,
        u_pe=u_pe,
        u_lw=u_lw,
        star_pos=star_pos,
    )


def radial_distance(pos, star_pos, boxsize):
    dx = pos - star_pos
    dx -= boxsize * np.round(dx / boxsize)
    return np.sqrt(np.sum(dx**2, axis=1))


def bin_radially(r, u, n_bins, r_max, empty=0.0):
    """Mean of `u` in `n_bins` equal radial bins spanning `0 .. r_max`.

    Particles outside the range are EXCLUDED, never clipped into the end
    bins: clipping mixes the far field into the outermost bin and biases
    its mean towards the far field's near-zero u.

    Parameters
    ----------
    r : numpy.ndarray
        Radius of each particle.
    u : numpy.ndarray
        Value to average, one per particle.
    n_bins : int
        Number of radial bins.
    r_max : float
        Outer edge of the binned range.
    empty : float, optional
        Value given to a bin holding no particle.

    Returns
    -------
    centres : numpy.ndarray
        Bin-centre radii.
    means : numpy.ndarray
        Per-bin mean of `u`, `empty` where the bin is unpopulated.
    """
    edges = np.linspace(0.0, r_max, n_bins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    idx = np.digitize(r, edges) - 1
    keep = (idx >= 0) & (idx < n_bins)
    idx, u_keep = idx[keep], u[keep]
    means = np.full(n_bins, empty)
    for i in range(n_bins):
        sel = idx == i
        if sel.sum() > 0:
            means[i] = u_keep[sel].mean()
    return centres, means


def outer_edge_above_threshold(r, u, threshold, n_bins, r_max):
    """Largest bin-centre radius where the radially-binned mean u still
    exceeds `threshold`; 0 if no bin exceeds it."""
    centres, means = bin_radially(r, u, n_bins, r_max, empty=0.0)
    above = means > threshold
    return float(centres[above].max()) if above.any() else 0.0, centres, means


def check_band(band, r, u, h_med, c_hyp, t, t0, opt, r_max):
    u_plateau = float(u.max())
    ok = True
    findings = []

    # Negativity: gated, unlike the Tier-1 script.
    epsilon_floor = 1e-6 * u_plateau if u_plateau > 0 else 1e-30
    n_neg = int(np.sum(u < -epsilon_floor))
    if n_neg > 0:
        ok = False
        findings.append(
            f"{band}: FAIL negativity -- {n_neg} particles below "
            f"-{epsilon_floor:.3e} (min u = {u.min():.3e})"
        )

    r_front = c_hyp * max(t - t0, 0.0)
    results = {}
    for eps in (0.1, 0.01, 0.001):
        r_edge, centres, means = outer_edge_above_threshold(
            r, u, eps * u_plateau, opt.n_bins, r_max
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

    # Bin-to-bin rise of u(r) behind the front, outside the near-source
    # exclusion zone. Reported only, never gated: see this module's
    # docstring. Bins below the smallest gated level carry no signal and
    # their ratio is noise over roughly zero, so they are left out.
    centres, means = bin_radially(r, u, opt.n_bins, r_max, empty=np.nan)
    exclude = centres < opt.near_source_h * h_med
    with np.errstate(invalid="ignore"):
        lit = means > 0.001 * u_plateau
    keep = ~np.isnan(means) & ~exclude & lit
    worst_bump = None
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
    check_finite("bulk time-step", np.atleast_1d(dt_bulk), opt.timesteps_log)
    print(f"Bulk gas time-step (timesteps.txt, n_gas={n_gas}): {dt_bulk:.6e}")

    all_ok = True
    fig, axes = plt.subplots(1, 2, figsize=(11, 5), sharey=True)
    colors = plt.cm.viridis(np.linspace(0, 0.9, len(files)))

    t0 = 0.0  # star BirthTime = 0 in this example's ICs (see README).

    for i, fn in enumerate(files):
        snap = load_snapshot(fn)
        t = snap["time"]
        if t <= t0:
            continue
        pos, h, ids = snap["pos"], snap["h"], snap["ids"]
        check_finite("Coordinates", pos, fn)
        check_finite("SmoothingLengths", h, fn)
        check_finite("star Coordinates", snap["star_pos"], fn)
        check_finite("BoxSize", np.atleast_1d(snap["boxsize"]), fn)
        check_finite("Time", np.atleast_1d(t), fn)
        gas_mask = np.ones(len(ids), dtype=bool)
        if opt.hot_particle_id >= 0:
            gas_mask = ids != opt.hot_particle_id
        h_med = float(np.median(h[gas_mask]))
        check_finite("median smoothing length h_med", np.atleast_1d(h_med), fn)
        c_hyp = min(opt.c_hyp_margin * h_med / dt_bulk, SPEED_OF_LIGHT_KM_S)
        check_finite("c_hyp", np.atleast_1d(c_hyp), fn)

        r_all = radial_distance(pos, snap["star_pos"], snap["boxsize"])
        # Full periodic minimum-image reach: no particle lies outside it, so
        # the gate can see a violation anywhere in the box.
        r_max = 0.5 * np.sqrt(3.0) * snap["boxsize"]

        fail_radius = c_hyp * max(t - t0, 0.0) + opt.fail_margin_h * h_med
        print(
            f"\n--- t={t:.4e}  h_med={h_med:.4e}  c_hyp={c_hyp:.4f}  "
            f"r_front={c_hyp * (t - t0):.4e} ({c_hyp * (t - t0) / h_med:.2f} h) ---"
        )
        if fail_radius >= r_max:
            all_ok = False
            print(
                f"  NOT EVALUABLE: the fail radius {fail_radius:.4e} is outside "
                f"the box's own reach {r_max:.4e}, so the causal-reach gate "
                f"cannot fire at this time. Shorten the run or enlarge the box."
            )

        for band, u_field in (("PE", "u_pe"), ("LW", "u_lw")):
            u_all = snap[u_field]
            check_finite(u_field, u_all, fn, band=band)
            r, u = r_all[gas_mask], u_all[gas_mask]
            res = check_band(band, r, u, h_med, c_hyp, t, t0, opt, r_max)
            all_ok &= res["ok"]
            status = "PASS" if res["ok"] else "FAIL"
            eps_str = ", ".join(
                f"eps={e}: edge={res['results'][e][0]/h_med:.2f}h "
                f"(C={res['results'][e][1]:.2f})"
                for e in (0.1, 0.01, 0.001)
            )
            bump = (
                "n/a (no lit bin outside the near-source zone)"
                if res["worst_bump"] is None
                else f"{res['worst_bump']:.3f} at r={res['worst_bump_r']:.4e}"
            )
            print(
                f"{band}: u_plateau={res['u_plateau']:.4e}  {eps_str}  "
                f"worst_bump (report only) = {bump}  -> {status}"
            )
            for finding in res["findings"]:
                print("  " + finding)

            ax = axes[0] if band == "PE" else axes[1]
            valid = np.isfinite(res["means"]) & (res["means"] > 0)
            ax.semilogy(
                res["centres"][valid] / h_med,
                res["means"][valid],
                "-",
                color=colors[i],
                label=f"t={t:.2e}",
            )
            ax.axvline(res["r_front"] / h_med, color=colors[i], ls="--", lw=1)

    for ax, band in zip(axes, ("PE", "LW")):
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
