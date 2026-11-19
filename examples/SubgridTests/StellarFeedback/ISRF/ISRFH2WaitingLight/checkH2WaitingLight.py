#!/usr/bin/env python3
"""Check the H2 line absorption of light that waits in a slot, in a clumpy cloud.

At a limited light speed the light that has not yet arrived waits in a reception slot
of the gas particle, and a partial release gives the released light the mean of the
whole pool (pro-rata rule). The light of one star in a static cloud reaches a particle
as shares from the neighbours that are closer to the star. Each share has the H2 column
of its sender plus the column of the hop. A pool releases a weighted mean of its shares,
so the released transmission and the released column lie between the smallest and the
largest value among the shares that the pool can hold. Two runs of one fixture, with the
light speed limited and with it unlimited (no waiting), mix the same shares with other
weights. Their values therefore differ by at most the ratio of that largest to that
smallest share value.

The reference is computed here and not by the code under test: the H2 density is
rebuilt from the snapshot particles as an SPH field on a grid, and the column of every
particle is the integral of that field along the straight line from the star to the
particle. The value of a share is the shielding function of the column of the sender
plus the mean H2 density of both ends of the hop times the hop length.

Checks, on the particles that both runs have lit and whose H2 density has changed by at
most 5 per cent since the first snapshot (the share values assume a frozen H2):
  hull_t   |ln(t_limited / t_control)| <= ln(largest / smallest share transmission)
  hull_N   |ln(N_limited / N_control)| <= ln(largest / smallest share column)
The control run must not wait. The limited run must lag the control run (its front is
behind) at the first compared snapshot, which shows that light waits in the limited run.

Exit codes: 0 passed, 1 a bar failed, 2 a file or field is missing, 3 a value is not
finite, 4 the run is void (a premise or a positive control is not met).
"""

import argparse
import os
import sys
from typing import List, Optional, Tuple

import numpy as np
from scipy.ndimage import map_coordinates
from scipy.spatial import cKDTree

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ISRFH2StarLit")
)
import checkH2Attenuation as G  # noqa: E402

GAMMA = G.KERNEL_GAMMA
GRID = 128
ROUND_OFF = 1e-5  # float32 rounding of the snapshot fields (relative)
H2_FROZEN = 0.05  # largest change of a particle's H2 density for the compared state
FIELDS = (
    "ParticleIDs",
    "Coordinates",
    "SmoothingLengths",
    "Densities",
    "InternalEnergies",
    "H2I",
    "Masses",
) + tuple(G.SPECIES_MASS_NUMBER)
LIT_FIELDS = FIELDS + G.REQUIRED_OWN + G.REQUIRED_COLUMN


def sph_field(snap: dict, weights: np.ndarray, n_grid: int) -> np.ndarray:
    """Scatter weights_j W(r, gamma h_j) of the particles onto a periodic grid of cell centres.

    W is the Wendland C2 kernel. The grid holds sum_j weights_j W.
    """
    box, pos, support = (
        snap["box"],
        snap["Coordinates"],
        GAMMA * snap["SmoothingLengths"],
    )
    dx = box / n_grid
    field = np.zeros((n_grid, n_grid, n_grid))
    ncell = int(np.ceil(support.max() / dx))
    index = np.floor(pos / dx).astype(int)
    span = range(-ncell, ncell + 1)
    for ox in span:
        for oy in span:
            for oz in span:
                cell = index + np.array([ox, oy, oz])
                q = np.linalg.norm((cell + 0.5) * dx - pos, axis=1) / support
                near = q < 1.0
                value = (
                    21.0
                    / (2.0 * np.pi * support[near] ** 3)
                    * (1.0 - q[near]) ** 4
                    * (1.0 + 4.0 * q[near])
                    * weights[near]
                )
                np.add.at(field, tuple((cell[near] % n_grid).T), value)
    return field


def line_integral(
    field: np.ndarray, star: np.ndarray, pos: np.ndarray, box: float
) -> Tuple[np.ndarray, np.ndarray]:
    """Return the integral of the grid field along the line from the star to each position.

    The line follows the nearest periodic image of the star. The field is sampled at
    half a cell along the line with linear interpolation. Units: field times box length.
    """
    dx = box / field.shape[0]
    vec = pos - star[None, :]
    vec -= box * np.round(vec / box)
    dist = np.linalg.norm(vec, axis=1)
    n_samples = int(np.ceil(dist.max() / (0.5 * dx))) + 1
    along = (np.arange(n_samples) + 0.5) / n_samples
    out = np.zeros(len(pos))
    for i in range(0, len(pos), 2048):
        sel = slice(i, i + 2048)
        pts = star[None, None, :] + vec[sel, None, :] * along[None, :, None]
        coord = (pts / dx - 0.5).reshape(-1, 3).T
        vals = map_coordinates(field, coord, order=1, mode="grid-wrap")
        out[sel] = vals.reshape(-1, n_samples).mean(axis=1) * dist[sel]
    return out, dist


def reference_columns(
    snap: dict, star: np.ndarray
) -> Tuple[np.ndarray, np.ndarray, float]:
    """Return the H2 column N_ref in cm^-2 of every particle, its distance and the field error.

    The field error is the median |field / n_H2 - 1| at the particle positions: the
    reconstructed SPH field must reproduce the particle values it is built from.
    """
    n_h2 = G.h2_column_density(snap)
    weights = snap["Masses"] / snap["Densities"] * n_h2
    field = sph_field(snap, weights, GRID)
    dx = snap["box"] / GRID
    at_particles = map_coordinates(
        field, (snap["Coordinates"] / dx - 0.5).T, order=1, mode="grid-wrap"
    )
    error = float(np.median(np.abs(at_particles / n_h2 - 1.0)))
    column, dist = line_integral(field, star, snap["Coordinates"], snap["box"])
    return column * snap["unit_length"], dist, error


def share_bounds(
    snap: dict,
    n_ref: np.ndarray,
    dist: np.ndarray,
    domain: np.ndarray,
    star: np.ndarray,
) -> Tuple[np.ndarray, ...]:
    """Return the smallest and largest share transmission and column of each domain particle.

    A share comes from a neighbour q closer to the star than the receiver p, within the
    support of q. It has the column N_ref(q) plus the mean H2 density of both ends times
    the hop length d_p - d_q. Its transmission is the shielding function at the receiver
    temperature and gas density. Particles without a closer neighbour get NaN.
    """
    box = snap["box"]
    pos = snap["Coordinates"]
    n_h2 = G.h2_column_density(snap)
    temp, dens = snap["Temperatures"], G.gas_number_density(snap)
    support = GAMMA * snap["SmoothingLengths"]
    tree = cKDTree(np.mod(pos, box), boxsize=box)
    idx = np.nonzero(domain)[0]
    lo_t = np.full(len(pos), np.nan)
    hi_t = np.full(len(pos), np.nan)
    lo_n = np.full(len(pos), np.nan)
    hi_n = np.full(len(pos), np.nan)
    for p, neigh in zip(
        idx, tree.query_ball_point(np.mod(pos[idx], box), support.max())
    ):
        neigh = np.asarray(neigh, dtype=int)
        d = pos[neigh] - pos[p]
        d -= box * np.round(d / box)
        up = (dist[neigh] < dist[p]) & (np.linalg.norm(d, axis=1) < support[neigh])
        if not up.any():
            continue
        q = neigh[up]
        hop = (dist[p] - dist[q]) * snap["unit_length"]
        col = n_ref[q] + 0.5 * (n_h2[p] + n_h2[q]) * hop
        t = G.transmission(col, np.full(len(q), temp[p]), np.full(len(q), dens[p]))
        lo_t[p], hi_t[p] = t.min(), t.max()
        lo_n[p], hi_n[p] = col.min(), col.max()
    return lo_t, hi_t, lo_n, hi_n


def lit_state(snap: dict) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the lit mask, the own transmission and the carried column in cm^-2."""
    own = snap["LWOwnSpecificEnergies"]
    dis = snap["LWOwnDissociatingSpecificEnergies"]
    lit = own > 0
    t = np.where(lit, dis / np.where(lit, own, 1.0), np.nan)
    return lit, t, snap[G.PARITY_COLUMN] / snap["unit_length"] ** 2


def load_run(run: str, number: int, first: dict) -> dict:
    snap = G.load_snapshot(
        os.path.join(run, "snap", f"snapshot_{number:04d}.hdf5"), LIT_FIELDS
    )
    G.align_to_ids(snap, first["ParticleIDs"], run)
    return snap


def run_state(snap: dict, star: np.ndarray, temps: np.ndarray):
    """Return the reference of one H2 state: column, distance, field error, domain and share bounds."""
    state = dict(snap)
    state["Temperatures"] = temps
    n_ref, dist, error = reference_columns(state, star)
    G.require_finite("reference", n_ref=n_ref, dist=dist)
    h_mean = float(snap["SmoothingLengths"].mean())
    domain = (dist > 2.0 * h_mean) & (dist < snap["box"] / 2.0 - 2.0 * h_mean)
    bounds = share_bounds(state, n_ref, dist, domain, star)
    return n_ref, dist, error, domain, bounds


def gate(args: argparse.Namespace) -> int:
    rep = G.Report()
    params = G.read_params(os.path.join(args.limited, "used_parameters.yml"))
    ctrl_params = G.read_params(os.path.join(args.control, "used_parameters.yml"))
    key = "GEARFeedback:ISRF_packet_light_speed"
    for name, p, ok in (
        ("limited run", params, ("fixed", "cap")),
        ("control run", ctrl_params, ("unlimited",)),
    ):
        got = str(p.get(key, "")).strip()
        if got not in ok:
            raise G.GateExit(
                G.EXIT_VOID, f"void run, the {name} has {key} = '{got}', not {ok}"
            )

    first_path = os.path.join(args.limited, "snap", "snapshot_0001.hdf5")
    first = G.load_snapshot(first_path, FIELDS)
    t_init = float(params.get("SPH:initial_temperature", "nan"))
    G.require_finite(
        "first snapshot", Temperatures=first["Temperatures"], t_init=np.array([t_init])
    )
    t_scale = t_init / float(np.median(first["Temperatures"]))
    star = G.star_particles(first_path)[0][0]

    ctl = load_run(args.control, args.control_step, first)
    lit_c, t_c, n_c = lit_state(ctl)
    ref_c = run_state(ctl, star, ctl["Temperatures"] * t_scale)
    domain = ref_c[3]
    saved = {"t_control": t_c, "N_control": n_c, "dist": ref_c[1], "domain": domain}
    rep.check(
        "reference field",
        ref_c[2] < 0.02,
        f"the rebuilt SPH H2 field reproduces the particle values to {ref_c[2]:.4f} median (bar 0.02)",
    )

    lag_checked = False
    for number in args.steps:
        lim = load_run(args.limited, number, first)
        lit_l, t_l, n_l = lit_state(lim)
        ref_l = run_state(lim, star, lim["Temperatures"] * t_scale)
        # The H2 only falls between the two snapshots, so the extreme share values of
        # the two states bound the values of every state in between.
        lo_t = np.fmin(ref_c[4][0], ref_l[4][0])
        hi_t = np.fmax(ref_c[4][1], ref_l[4][1])
        lo_n = np.fmin(ref_c[4][2], ref_l[4][2])
        hi_n = np.fmax(ref_c[4][3], ref_l[4][3])
        with_shares = domain & np.isfinite(lo_t)
        G.require_finite("share bounds", lo=lo_t[with_shares], hi=hi_t[with_shares])
        width_t = np.log(hi_t / lo_t)
        width_n = np.log(hi_n / lo_n)
        n_first = G.h2_column_density(first)
        frozen = (np.abs(G.h2_column_density(lim) / n_first - 1.0) <= H2_FROZEN) & (
            np.abs(G.h2_column_density(ctl) / n_first - 1.0) <= H2_FROZEN
        )
        both = lit_l & lit_c & with_shares
        n_before = int(both.sum())
        both &= frozen
        if both.sum() < 0.5 * n_before:
            raise G.GateExit(
                G.EXIT_VOID,
                f"void run, at snapshot {number} the H2 of {n_before - int(both.sum())} of "
                f"{n_before} particles changed by more than {H2_FROZEN:g}: the H2 is not frozen",
            )
        rep.info(
            f"frozen H2, step {number}",
            f"{n_before - int(both.sum())} of {n_before} particles left out (H2 changed by more than {H2_FROZEN:g})",
        )
        saved[f"t_limited_{number}"], saved[f"N_limited_{number}"] = t_l, n_l
        saved[f"width_t_{number}"], saved[f"width_N_{number}"] = width_t, width_n
        saved[f"t_ref_{number}"] = G.transmission(
            ref_l[0], lim["Temperatures"] * t_scale, G.gas_number_density(lim)
        )
        if not lag_checked:
            lag_checked = True
            n_lim, n_ctl = int((lit_l & domain).sum()), int((lit_c & domain).sum())
            if not n_lim < n_ctl:
                raise G.GateExit(
                    G.EXIT_VOID,
                    f"void run, at snapshot {number} the limited run lit {n_lim} domain "
                    f"particles and the control {n_ctl}: no light waits in the limited run",
                )
            rep.info(
                "positive control",
                f"snapshot {number}: limited run lit {n_lim}, control {n_ctl}",
            )
        G.require_finite(
            f"snapshot {number}",
            t_limited=t_l[both],
            t_control=t_c[both],
            N_limited=n_l[both],
            N_control=n_c[both],
        )
        if not both.any() or not np.all(t_l[both] > 0) or not np.all(t_c[both] > 0):
            raise G.GateExit(
                G.EXIT_VOID, f"void run, snapshot {number} has no usable lit particle"
            )
        d_t = np.abs(np.log(t_l[both] / t_c[both]))
        positive = (n_l[both] > 0) & (n_c[both] > 0)
        d_n = np.full(int(both.sum()), np.inf)
        d_n[positive] = np.abs(np.log(n_l[both][positive] / n_c[both][positive]))
        over_t = d_t - width_t[both] - ROUND_OFF
        over_n = d_n - width_n[both] - ROUND_OFF
        t_ref = saved[f"t_ref_{number}"]
        e_l, e_c = np.log(t_l[both] / t_ref[both]), np.log(t_c[both] / t_ref[both])
        rep.info(
            f"error against the ray, step {number}",
            f"ln(t/t_ray) median [10, 90 per cent]: limited {np.median(e_l):+.3f} "
            f"[{np.percentile(e_l, 10):+.3f}, {np.percentile(e_l, 90):+.3f}], control "
            f"{np.median(e_c):+.3f} [{np.percentile(e_c, 10):+.3f}, {np.percentile(e_c, 90):+.3f}]",
        )
        rep.check(
            f"hull_t step {number}",
            bool(np.all(over_t <= 0)),
            f"{int(both.sum())} particles; largest |ln(t_lim/t_ctl)| = {d_t.max():.4f}; "
            f"allowed spread median {np.median(width_t[both]):.3f}; "
            f"median fraction of it {np.median(d_t / np.maximum(width_t[both], 1e-12)):.3f}; "
            f"worst excess over the bar {over_t.max():.2e}",
        )
        rep.check(
            f"hull_N step {number}",
            bool(np.all(over_n <= 0)),
            f"largest |ln(N_lim/N_ctl)| = {d_n.max():.4f}; worst excess over the bar {over_n.max():.2e}",
        )
    if args.save:
        np.savez(
            args.save,
            steps=np.array(args.steps),
            control_step=args.control_step,
            **saved,
        )
    return rep.exit_code()


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    p.add_argument(
        "--limited", required=True, help="Run directory with the limited light speed"
    )
    p.add_argument(
        "--control", required=True, help="Run directory with the unlimited light speed"
    )
    p.add_argument(
        "--control-step",
        required=True,
        type=int,
        help="Snapshot of the control run: the front has crossed the domain and the "
        "transport has settled, before the H2 changes",
    )
    p.add_argument(
        "--steps",
        required=True,
        type=lambda s: [int(x) for x in s.split(",")],
        help="Comma-separated snapshot numbers of the limited run to compare",
    )
    p.add_argument(
        "--save", help="Write the arrays of the comparison to this .npz file"
    )
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        code = gate(args)
    except G.GateExit as stop:
        print(f"FAIL  {stop.message}")
        code = stop.code
    print(f"exit class {code}")
    return code


if __name__ == "__main__":
    sys.exit(main())
