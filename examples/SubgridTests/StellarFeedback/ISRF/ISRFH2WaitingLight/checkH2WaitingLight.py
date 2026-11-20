#!/usr/bin/env python3
"""Check the H2 line absorption of light that waits in a slot, in a clumpy cloud.

At a limited light speed the light that has not yet arrived waits in a reception slot
of the gas particle, and a partial release gives the released light the mean of the
whole pool (pro-rata rule). The light of one star in a static cloud reaches a particle
as shares from the neighbours that are closer to the star. Each share has the H2 column
carried by its sender plus the column of the hop. A pool releases a weighted mean of its
shares, so the transmission and the column that a particle carries lie between the
smallest and the largest value among the shares that its pool held.

Two runs of one fixture, with the light speed limited and unlimited (no waiting), are
read at the same snapshot time. Each run has its own shares, built from the carried
columns of its own senders in every snapshot up to the compared one, because light can
wait in a slot for an unknown number of steps before it is released.
Both runs mix shares of their own set, so their values differ by at most
ln(largest / smallest value over the shares of both runs).

Checks, on the particles that both runs have lit and whose H2 density has changed by at
most 5 per cent since the first snapshot:
  hull_t   |ln(t_limited / t_control)| <= ln(largest / smallest share transmission)
  hull_N   min share column <= N_limited, N_control <= max share column
  own_t, own_N   the limited run alone lies inside the range of its own shares
The control run must not wait. The limited run must lag the control run (its front is
behind) at the first compared snapshot, which shows that light waits in the limited run.

No snapshot index is chosen by the user: every snapshot that both runs hold is read, in
time order, from the first at which at least MIN_COMPARE particles are lit in both. The
control run must itself obey the range of its own shares (the premise of the shares); the
comparison ends at the first snapshot where it does not, because the carried column of
the control run then no longer follows the shares. The comparison must reach at least one
snapshot, else the run is void.

Exit codes: 0 passed, 1 a bar failed, 2 a file or field is missing, 3 a value is not
finite, 4 the run is void (a premise or a positive control is not met).
"""

import argparse
import glob
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ISRFH2StarLit")
)
import checkH2Attenuation as G  # noqa: E402

GAMMA = G.KERNEL_GAMMA
ROUND_OFF = 1e-5  # float32 rounding of the snapshot fields (relative)
H2_FROZEN = 0.05  # largest change of a particle's H2 density for the compared state
MIN_COMPARE = 200  # fewest particles lit in both runs for a snapshot to be compared
TIME_MATCH = 1e-9  # relative tolerance on the snapshot times of the two runs
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


def lit_state(snap: dict) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the lit mask, the own transmission and the carried column in cm^-2."""
    own = snap["LWOwnSpecificEnergies"]
    dis = snap["LWOwnDissociatingSpecificEnergies"]
    lit = own > 0
    t = np.where(lit, dis / np.where(lit, own, 1.0), np.nan)
    return lit, t, snap[G.PARITY_COLUMN] / snap["unit_length"] ** 2


class Run:
    """Snapshots of one run, read in the particle order of a reference snapshot."""

    def __init__(self, path: str, first: dict, t_scale: float) -> None:
        self.path, self.first, self.t_scale = path, first, t_scale
        self.cache: Dict[int, dict] = {}

    def get(self, number: int) -> dict:
        """Return snapshot number with the lit state, the temperature and the H2 density."""
        if number not in self.cache:
            snap = G.load_snapshot(
                os.path.join(self.path, "snap", f"snapshot_{number:04d}.hdf5"),
                LIT_FIELDS,
            )
            G.align_to_ids(snap, self.first["ParticleIDs"], self.path)
            snap["Temperatures"] = snap["Temperatures"] * self.t_scale
            snap["lit"], snap["t"], snap["N"] = lit_state(snap)
            snap["n_h2"] = G.h2_column_density(snap)
            snap["n_gas"] = G.gas_number_density(snap)
            self.cache[number] = snap
        return self.cache[number]


def upstream_neighbours(
    snap: dict, star: np.ndarray, domain: np.ndarray, dist: np.ndarray
) -> List[Tuple[int, np.ndarray]]:
    """Return, for each domain particle p, the neighbours q closer to the star than p.

    A neighbour q sends a share to p when p lies within the support of q.
    """
    box, pos = snap["box"], snap["Coordinates"]
    support = GAMMA * snap["SmoothingLengths"]
    tree = cKDTree(np.mod(pos, box), boxsize=box)
    idx = np.nonzero(domain)[0]
    out = []
    for p, neigh in zip(
        idx, tree.query_ball_point(np.mod(pos[idx], box), support.max())
    ):
        neigh = np.asarray(neigh, dtype=int)
        d = pos[neigh] - pos[p]
        d -= box * np.round(d / box)
        up = (dist[neigh] < dist[p]) & (np.linalg.norm(d, axis=1) < support[neigh])
        out.append((p, neigh[up]))
    return out


def share_bounds(
    run: Run,
    number: int,
    neighbours: List[Tuple[int, np.ndarray]],
    dist: np.ndarray,
) -> Tuple[np.ndarray, ...]:
    """Return the smallest and largest share transmission and column of each particle.

    A share from q has the column carried by q plus the mean H2 density of both ends times
    the hop length d_p - d_q, taken from every snapshot up to the compared one at which q
    is lit (light may have waited since then). Its transmission is the shielding function at the receiver temperature and
    gas density. Particles without a lit upstream neighbour get NaN.
    """
    snap = run.get(number)
    window = [run.get(j) for j in range(1, number + 1)]
    size = len(snap["ParticleIDs"])
    lo_t, hi_t, lo_n, hi_n = (np.full(size, np.nan) for _ in range(4))
    for p, q in neighbours:
        if len(q) == 0:
            continue
        hop = (dist[p] - dist[q]) * snap["unit_length"]
        cols = []
        for past in window:
            sent = past["lit"][q]
            if sent.any():
                cols.append(
                    past["N"][q[sent]]
                    + 0.5 * (snap["n_h2"][p] + snap["n_h2"][q[sent]]) * hop[sent]
                )
        if not cols:
            continue
        col = np.concatenate(cols)
        t = G.transmission(
            col,
            np.full(len(col), snap["Temperatures"][p]),
            np.full(len(col), snap["n_gas"][p]),
        )
        lo_t[p], hi_t[p], lo_n[p], hi_n[p] = t.min(), t.max(), col.min(), col.max()
    return lo_t, hi_t, lo_n, hi_n


def outside(value: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Return the mask of values outside [lo, hi] by more than the float rounding."""
    return (value > hi * (1.0 + ROUND_OFF)) | (value < lo * (1.0 - ROUND_OFF))


def snapshot_numbers(run: str) -> List[int]:
    """Return the snapshot numbers found in the snap directory of a run, sorted."""
    names = glob.glob(os.path.join(run, "snap", "snapshot_*.hdf5"))
    return sorted(
        int(m.group(1)) for n in names if (m := re.search(r"_(\d+)\.hdf5$", n))
    )


def check_light_speed(args: argparse.Namespace) -> dict:
    """Return the parameters of the limited run; stop if either run has the wrong light speed."""
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
    return params


def gate(args: argparse.Namespace) -> int:
    """Compare the limited and the unlimited run at equal snapshot times; return the exit class."""
    rep = G.Report()
    params = check_light_speed(args)
    first_path = os.path.join(args.limited, "snap", "snapshot_0001.hdf5")
    first = G.load_snapshot(first_path, FIELDS)
    t_init = float(params.get("SPH:initial_temperature", "nan"))
    G.require_finite(
        "first snapshot", Temperatures=first["Temperatures"], t_init=np.array([t_init])
    )
    t_scale = t_init / float(np.median(first["Temperatures"]))
    star = G.star_particles(first_path)[0][0]
    lim_run, ctl_run = Run(args.limited, first, t_scale), Run(
        args.control, first, t_scale
    )

    box = first["box"]
    vec = first["Coordinates"] - star[None, :]
    vec -= box * np.round(vec / box)
    dist = np.linalg.norm(vec, axis=1)
    h_mean = float(first["SmoothingLengths"].mean())
    domain = (dist > 2.0 * h_mean) & (dist < box / 2.0 - 2.0 * h_mean)
    n_first = G.h2_column_density(first)

    numbers = sorted(
        set(snapshot_numbers(args.limited)) & set(snapshot_numbers(args.control))
    )
    compared: List[int] = []
    saved: Dict[str, np.ndarray] = {}
    ended = None
    for number in [n for n in numbers if n >= 2]:
        lim, ctl = lim_run.get(number), ctl_run.get(number)
        if not np.isclose(lim["time"], ctl["time"], rtol=TIME_MATCH, atol=0.0):
            raise G.GateExit(
                G.EXIT_VOID,
                f"void run, snapshot {number} has the time {lim['time']:.9e} in the limited "
                f"run and {ctl['time']:.9e} in the control run",
            )
        both = domain & lim["lit"] & ctl["lit"]
        if int(both.sum()) < MIN_COMPARE:
            continue
        neighbours = upstream_neighbours(lim, star, domain, dist)
        b_lim = share_bounds(lim_run, number, neighbours, dist)
        b_ctl = share_bounds(ctl_run, number, neighbours, dist)
        with_shares = np.isfinite(b_lim[0]) & np.isfinite(b_ctl[0])
        G.require_finite(
            "share bounds",
            lo=np.minimum(b_lim[0], b_ctl[0])[both & with_shares],
            hi=np.maximum(b_lim[1], b_ctl[1])[both & with_shares],
        )
        steady = both & np.isfinite(b_ctl[0])
        broken = outside(
            ctl["t"][steady], b_ctl[0][steady], b_ctl[1][steady]
        ) | outside(ctl["N"][steady], b_ctl[2][steady], b_ctl[3][steady])
        if broken.any():
            ended = (number, int(broken.sum()), int(steady.sum()))
            break
        lo_t, hi_t = np.fmin(b_lim[0], b_ctl[0]), np.fmax(b_lim[1], b_ctl[1])
        lo_n, hi_n = np.fmin(b_lim[2], b_ctl[2]), np.fmax(b_lim[3], b_ctl[3])
        frozen = (np.abs(lim["n_h2"] / n_first - 1.0) <= H2_FROZEN) & (
            np.abs(ctl["n_h2"] / n_first - 1.0) <= H2_FROZEN
        )
        use = both & with_shares
        n_before = int(use.sum())
        use &= frozen
        if use.sum() < 0.5 * n_before:
            raise G.GateExit(
                G.EXIT_VOID,
                f"void run, at snapshot {number} the H2 of {n_before - int(use.sum())} of "
                f"{n_before} particles changed by more than {H2_FROZEN:g}: the H2 is not frozen",
            )
        if not compared:
            n_lim, n_ctl = int((lim["lit"] & domain).sum()), int(
                (ctl["lit"] & domain).sum()
            )
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
        compared.append(number)
        G.require_finite(
            f"snapshot {number}",
            t_limited=lim["t"][use],
            t_control=ctl["t"][use],
            N_limited=lim["N"][use],
            N_control=ctl["N"][use],
        )
        if not np.all(lim["t"][use] > 0) or not np.all(ctl["t"][use] > 0):
            raise G.GateExit(
                G.EXIT_VOID,
                f"void run, snapshot {number} has a lit particle with t = 0",
            )
        width_t = np.log(hi_t[use] / lo_t[use])
        d_t = np.abs(np.log(lim["t"][use] / ctl["t"][use]))
        over_t = d_t - width_t - ROUND_OFF
        n_out = outside(lim["N"][use], lo_n[use], hi_n[use]) | outside(
            ctl["N"][use], lo_n[use], hi_n[use]
        )
        own_t = outside(lim["t"][use], b_lim[0][use], b_lim[1][use])
        own_n = outside(lim["N"][use], b_lim[2][use], b_lim[3][use])
        rep.check(
            f"hull_t snapshot {number}",
            bool(np.all(over_t <= 0)),
            f"{int(use.sum())} particles; largest |ln(t_lim/t_ctl)| = {d_t.max():.4f}; "
            f"allowed spread min {width_t.min():.3f}, median {np.median(width_t):.3f}; "
            f"worst excess over the bar {over_t.max():.2e}",
        )
        rep.check(
            f"hull_N snapshot {number}",
            not n_out.any(),
            f"{int(n_out.sum())} of {int(use.sum())} particles carry a column outside the range of the shares",
        )
        rep.check(
            f"own range snapshot {number}",
            not (own_t.any() or own_n.any()),
            f"limited run: {int(own_t.sum())} transmissions, {int(own_n.sum())} columns outside its own shares",
        )
        saved[f"d_t_{number}"] = d_t
        saved[f"width_t_{number}"] = width_t
        saved[f"own_excess_t_{number}"] = np.maximum(
            np.log(lim["t"][use] / b_lim[1][use]), np.log(b_lim[0][use] / lim["t"][use])
        )
    if ended is not None:
        rep.info(
            "end of comparison",
            f"at snapshot {ended[0]} the control run has {ended[1]} of {ended[2]} lit particles "
            "outside the range of its own shares: later snapshots are not compared",
        )
    if not compared:
        raise G.GateExit(
            G.EXIT_VOID,
            f"void run, no snapshot has {MIN_COMPARE} particles lit in both runs with a steady control run",
        )
    rep.info("snapshots compared", ", ".join(str(n) for n in compared))
    if args.save:
        np.savez(args.save, snapshots=np.array(compared), **saved)
    return rep.exit_code()


def build_parser() -> argparse.ArgumentParser:
    """Return the command-line parser."""
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
        "--save", help="Write the arrays of the comparison to this .npz file"
    )
    return p


def main(argv: Optional[List[str]] = None) -> int:
    """Run the check and return the exit class."""
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
