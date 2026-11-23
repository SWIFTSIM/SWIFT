#!/usr/bin/env python3
"""Check the H2 line absorption of light that waits in a slot, in a clumpy cloud.

At a limited light speed the light that has not yet arrived waits in a reception slot
of the gas particle, and a partial release gives the released light the mean of the
whole pool (pro-rata rule). A pool releases a weighted mean of its shares, so the
transmission and the column that a particle carries lie between the smallest and the
largest value among the shares that its pool held.

Definitions. The light of one star in a static cloud reaches a particle as shares from
the neighbours that are closer to the star (the senders). The share of a sender has the
H2 column N carried by the sender plus the column dN of the hop between the two, and the
transmission carried by the sender times f(N + dN) / f(N), with f the shielding function
at the temperature and density of the receiver. A sender value of snapshot j was sent
between snapshots j - 1 and j: the hop is built with the H2 of both snapshots and both
shares enter the range. The range of a particle is the smallest to the largest share
value. The shares come from the carried values of the senders of the same run in all
snapshots up to the compared one, because light can wait in a slot for an unknown number
of steps before it is released.

Two runs of one fixture, with the light speed limited and unlimited (no waiting), are
read at the same snapshot time. Both runs mix shares of their own set, so their values
differ by at most ln(largest / smallest value over the shares of both runs).

Checks, on the particles that both runs have lit and whose H2 density has changed by at
most 5 per cent since the first snapshot:
  hull_t   |ln(t_limited / t_control)| <= ln(largest / smallest share transmission)
  hull_N   min share column <= N_limited, N_control <= max share column
  own_t, own_N   the limited run alone lies inside the range of its own shares
and, on all lit domain particles of the control run, without the H2 condition:
  control        the control run lies inside the range of its own shares
The control run must not wait. The limited run must lag the control run (its front is
behind) at the first compared snapshot, which shows that light waits in the limited run.

No snapshot index is chosen by the user: every snapshot that both runs hold is read, in
time order, from the first at which at least MIN_COMPARE particles are lit in both. The
control test starts at that first compared snapshot. Both runs must be non-periodic: the
light of the periodic images of the star reaches a particle through neighbours that are
not shares. The control run has no waiting light and shares the transport code with the
limited run, so on all its lit domain particles it must lie inside the range of its own
shares. If it does not, the check fails (exit 1) at the first such snapshot and reads no
later one. The comparison ends without a failure where the H2 of more than half the
compared particles has changed. The check prints the first and last snapshot read and why
the comparison ends. The comparison must reach at least one snapshot, else the run is
void.

Finiteness: when a snapshot is read, the transmission and the column of every lit
particle (inside the domain or not), the H2 density, the gas density and the temperature
of every particle are tested, before any share is built. A value that is not finite stops
the check (exit 3). Particles without a lit upstream neighbour are marked by a mask.

Exit codes: 0 passed, 1 a bar failed, 2 a file or field is missing, 3 a value is not
finite, 4 the run is void (a premise or a positive control is not met).

Scope and limits.
A pass means that, at every snapshot read, each compared value lies within the smallest to
largest range of the shares that the check rebuilds from the carried values of its
upstream neighbours, at any snapshot up to the compared one. It does not see an error that
is the same on every particle, in either sign. It does not see a pooling error that keeps
a weighted mean of the shares: a lost share, a share counted twice, wrong weights, or
light released twice. The comparison between the two runs follows from the two own-range
tests and is not separate evidence. The check does not compare a column with the straight
line from the star. It sees an error that pushes a particle outside the range, for example
an H2 shielding that is too strong or too weak on each hop. An error of the transport code
that both runs share is seen only through the control run, and only when it leaves the
range of its own shares.
The control run can also leave its range for a reason of the fixture: a sender that is lit
between two snapshots at the light front is not in the shares. A change of `level` or
`gas_density` can move such a front into the compared window, and the check then fails
without a defect of the code.
The `cap` mode is refused (exit 4). The snapshots hold the mean of the slot columns of a
particle, not the column of the slot that a share leaves from, no slot origin and no send
time. A value that sits on the edge of its range because it has a single share can then
lie outside it, and the check cannot tell this from a defect.
"""

import argparse
import glob
import os
import re
import sys
from typing import Dict, List, NamedTuple, Optional, Tuple

import numpy as np
from scipy.spatial import cKDTree

sys.path.insert(
    0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "ISRFH2StarLit")
)
import checkH2Attenuation as attenuation  # noqa: E402

GAMMA = attenuation.KERNEL_GAMMA
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
) + tuple(attenuation.SPECIES_MASS_NUMBER)
LIT_FIELDS = FIELDS + attenuation.REQUIRED_OWN + attenuation.REQUIRED_COLUMN


class ShareRange(NamedTuple):
    """Smallest and largest share transmission and column of each particle."""

    t_min: np.ndarray
    t_max: np.ndarray
    n_min: np.ndarray
    n_max: np.ndarray
    has_shares: np.ndarray  # False where no lit upstream neighbour sends a share


class Cloud(NamedTuple):
    """Geometry of the fixture, from the first snapshot."""

    star: np.ndarray
    dist: np.ndarray  # distance of each particle to the star
    domain: np.ndarray  # mask of the particles that the check reads
    n_h2_first: np.ndarray  # H2 density of each particle at the first snapshot


def lit_state(snap: dict) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the lit mask, the own transmission and the carried column in cm^-2."""
    own = snap["LWOwnSpecificEnergies"]
    dis = snap["LWOwnDissociatingSpecificEnergies"]
    lit = own > 0
    t = np.where(lit, dis / np.where(lit, own, 1.0), np.nan)
    return lit, t, snap[attenuation.PARITY_COLUMN] / snap["unit_length"] ** 2


class Run:
    """Snapshots of one run, read in the particle order of a reference snapshot."""

    def __init__(self, path: str, first: dict, t_scale: float) -> None:
        self.path, self.first, self.t_scale = path, first, t_scale
        self.cache: Dict[int, dict] = {}

    def get(self, number: int) -> dict:
        """Return snapshot number with the lit state, the temperature and the H2 density.

        Stop the check when a value that a share uses is not finite.
        """
        if number not in self.cache:
            snap = attenuation.load_snapshot(
                os.path.join(self.path, "snap", f"snapshot_{number:04d}.hdf5"),
                LIT_FIELDS,
            )
            attenuation.align_to_ids(snap, self.first["ParticleIDs"], self.path)
            snap["Temperatures"] = snap["Temperatures"] * self.t_scale
            snap["lit"], snap["t"], snap["N"] = lit_state(snap)
            snap["n_h2"] = attenuation.h2_column_density(snap)
            snap["n_gas"] = attenuation.gas_number_density(snap)
            lit = snap["lit"]
            attenuation.require_finite(
                f"{self.path}, snapshot {number}",
                t=snap["t"][lit],
                N=snap["N"][lit],
                n_h2=snap["n_h2"],
                n_gas=snap["n_gas"],
                temperature=snap["Temperatures"],
            )
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
) -> ShareRange:
    """Return the range of the shares that each particle receives, at snapshot number.

    The column of a share is the column carried by the sender q plus the mean H2 density
    of both ends times the hop length d_p - d_q. The transmission is the one carried by
    q times f_p(N_q + dN) / f_p(N_q), with f_p the shielding function at the receiver
    temperature and gas density. The values of q come from every snapshot up to the
    compared one at which q is lit (light may have waited since then). The H2 of the hop
    is read at j - 1 and at j for the sender value of snapshot j.
    """
    snap = run.get(number)
    size = len(snap["ParticleIDs"])
    t_min, t_max, n_min, n_max = (np.full(size, np.nan) for _ in range(4))
    has_shares = np.zeros(size, dtype=bool)
    for p, q in neighbours:
        if len(q) == 0:
            continue
        hop = (dist[p] - dist[q]) * snap["unit_length"]
        cols, trans = [], []
        for j in range(1, number + 1):
            past = run.get(j)
            sent = past["lit"][q]
            if not sent.any():
                continue
            senders = q[sent]
            col_sender = past["N"][senders]
            temp = np.full(len(senders), snap["Temperatures"][p])
            dens = np.full(len(senders), snap["n_gas"][p])
            t_sender = past["t"][senders]
            f_sender = attenuation.transmission(col_sender, temp, dens)
            for hop_snap in [past] if j == 1 else [run.get(j - 1), past]:
                col_hop = (
                    0.5 * (hop_snap["n_h2"][p] + hop_snap["n_h2"][senders]) * hop[sent]
                )
                cols.append(col_sender + col_hop)
                trans.append(
                    t_sender
                    * attenuation.transmission(col_sender + col_hop, temp, dens)
                    / f_sender
                )
        if not cols:
            continue
        col, t = np.concatenate(cols), np.concatenate(trans)
        t_min[p], t_max[p], n_min[p], n_max[p] = t.min(), t.max(), col.min(), col.max()
        has_shares[p] = True
    return ShareRange(t_min, t_max, n_min, n_max, has_shares)


def outside(value: np.ndarray, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Return the mask of values outside [lo, hi] by more than the float rounding."""
    return (value > hi * (1.0 + ROUND_OFF)) | (value < lo * (1.0 - ROUND_OFF))


def snapshot_numbers(run: str) -> List[int]:
    """Return the snapshot numbers found in the snap directory of a run, sorted."""
    names = glob.glob(os.path.join(run, "snap", "snapshot_*.hdf5"))
    return sorted(
        int(m.group(1)) for n in names if (m := re.search(r"_(\d+)\.hdf5$", n))
    )


def check_premises(args: argparse.Namespace) -> dict:
    """Return the parameters of the limited run; stop if a run does not meet a premise.

    The limited run must use the fixed light speed, the control run the unlimited one,
    and both runs must be non-periodic.
    """
    params = attenuation.read_params(os.path.join(args.limited, "used_parameters.yml"))
    ctrl_params = attenuation.read_params(
        os.path.join(args.control, "used_parameters.yml")
    )
    key = "GEARFeedback:ISRF_packet_light_speed"
    for name, p, ok in (
        ("limited run", params, ("fixed",)),
        ("control run", ctrl_params, ("unlimited",)),
    ):
        got = str(p.get(key, "")).strip()
        if got == "cap":
            raise attenuation.GateExit(
                attenuation.EXIT_VOID,
                f"void run, the {name} uses the cap mode, which this check does not support: "
                "the snapshots hold neither the slot column, the slot origin nor the send time "
                "of a share, so the check cannot bound them",
            )
        if got not in ok:
            raise attenuation.GateExit(
                attenuation.EXIT_VOID,
                f"void run, the {name} has {key} = '{got}', not {ok}",
            )
        periodic = str(p.get("InitialConditions:periodic", "")).strip()
        if periodic != "0":
            raise attenuation.GateExit(
                attenuation.EXIT_VOID,
                f"void run, the {name} has InitialConditions:periodic = '{periodic}': "
                "light of the periodic images of the star is not in the shares",
            )
    return params


def locate_cloud(first: dict, star: np.ndarray) -> Cloud:
    """Return the distances to the star and the domain that the check reads.

    The domain runs from 2 smoothing lengths to half a box minus 2 smoothing lengths.
    """
    box = first["box"]
    vec = first["Coordinates"] - star[None, :]
    vec -= box * np.round(vec / box)
    dist = np.linalg.norm(vec, axis=1)
    h_mean = float(first["SmoothingLengths"].mean())
    domain = (dist > 2.0 * h_mean) & (dist < box / 2.0 - 2.0 * h_mean)
    return Cloud(star, dist, domain, attenuation.h2_column_density(first))


def control_leaves_range(
    rep: attenuation.Report,
    number: int,
    ctl: dict,
    shares: ShareRange,
    cloud: Cloud,
) -> Optional[str]:
    """Fail the check if the control run is outside the range of its own shares.

    Return the message of the failure, or None when the control run is inside.
    """
    checked = cloud.domain & ctl["lit"] & shares.has_shares
    attenuation.require_finite(
        f"control run, snapshot {number}, share bounds",
        t_min=shares.t_min[checked],
        t_max=shares.t_max[checked],
        n_min=shares.n_min[checked],
        n_max=shares.n_max[checked],
    )
    miss = outside(
        ctl["t"][checked], shares.t_min[checked], shares.t_max[checked]
    ) | outside(ctl["N"][checked], shares.n_min[checked], shares.n_max[checked])
    if not miss.any():
        return None
    message = f"the control run leaves the range of its own shares at snapshot {number}"
    rep.check(
        f"control range snapshot {number}",
        False,
        f"{message}: {int(miss.sum())} of {int(checked.sum())} lit domain particles are "
        "outside; the control has no waiting light, so the transport code that both runs "
        "share does not follow the hop rule, or the fixture moves a light front into the window",
    )
    return message


def compare_snapshot(
    rep: attenuation.Report,
    number: int,
    runs: Tuple[dict, dict],
    shares: Tuple[ShareRange, ShareRange],
    both: np.ndarray,
    cloud: Cloud,
    first_compared: bool,
    saved: Dict[str, np.ndarray],
) -> Optional[str]:
    """Test the limited run and the two runs together at one snapshot.

    Return why the comparison ends at this snapshot (the H2 has changed), else None.
    Stop the check when the run is void.
    """
    lim, ctl = runs
    lim_shares, ctl_shares = shares
    use = both & lim_shares.has_shares & ctl_shares.has_shares
    lo_t = np.fmin(lim_shares.t_min, ctl_shares.t_min)
    hi_t = np.fmax(lim_shares.t_max, ctl_shares.t_max)
    lo_n = np.fmin(lim_shares.n_min, ctl_shares.n_min)
    hi_n = np.fmax(lim_shares.n_max, ctl_shares.n_max)
    attenuation.require_finite(
        "share bounds", lo=lo_t[use], hi=hi_t[use], lo_n=lo_n[use], hi_n=hi_n[use]
    )
    frozen = (np.abs(lim["n_h2"] / cloud.n_h2_first - 1.0) <= H2_FROZEN) & (
        np.abs(ctl["n_h2"] / cloud.n_h2_first - 1.0) <= H2_FROZEN
    )
    n_before = int(use.sum())
    use &= frozen
    changed = n_before - int(use.sum())
    if int(use.sum()) < 0.5 * n_before:
        if first_compared:
            raise attenuation.GateExit(
                attenuation.EXIT_VOID,
                f"void run, at snapshot {number} the H2 of {changed} of {n_before} "
                f"particles changed by more than {H2_FROZEN:g}: the H2 is not frozen",
            )
        return (
            f"at snapshot {number} the H2 of {changed} of {n_before} particles changed "
            f"by more than {H2_FROZEN:g}: this and later snapshots are not read"
        )
    if first_compared:
        n_lim = int((lim["lit"] & cloud.domain).sum())
        n_ctl = int((ctl["lit"] & cloud.domain).sum())
        if not n_lim < n_ctl:
            raise attenuation.GateExit(
                attenuation.EXIT_VOID,
                f"void run, at snapshot {number} the limited run lit {n_lim} domain "
                f"particles and the control {n_ctl}: no light waits in the limited run",
            )
        rep.info(
            "positive control",
            f"snapshot {number}: limited run lit {n_lim}, control {n_ctl}",
        )
    attenuation.require_finite(
        f"snapshot {number}",
        t_limited=lim["t"][use],
        t_control=ctl["t"][use],
        N_limited=lim["N"][use],
        N_control=ctl["N"][use],
    )
    if not np.all(lim["t"][use] > 0) or not np.all(ctl["t"][use] > 0):
        raise attenuation.GateExit(
            attenuation.EXIT_VOID,
            f"void run, snapshot {number} has a lit particle with t = 0",
        )
    width_t = np.log(hi_t[use] / lo_t[use])
    d_t = np.abs(np.log(lim["t"][use] / ctl["t"][use]))
    over_t = d_t - width_t - ROUND_OFF
    n_out = outside(lim["N"][use], lo_n[use], hi_n[use]) | outside(
        ctl["N"][use], lo_n[use], hi_n[use]
    )
    own_t = outside(lim["t"][use], lim_shares.t_min[use], lim_shares.t_max[use])
    own_n = outside(lim["N"][use], lim_shares.n_min[use], lim_shares.n_max[use])
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
        f"{int(n_out.sum())} of {int(use.sum())} particles carry a column outside "
        "the range of the shares",
    )
    rep.check(
        f"own range snapshot {number}",
        not (own_t.any() or own_n.any()),
        f"limited run: {int(own_t.sum())} transmissions, {int(own_n.sum())} columns "
        "outside its own shares",
    )
    saved[f"d_t_{number}"] = d_t
    saved[f"width_t_{number}"] = width_t
    saved[f"own_excess_t_{number}"] = np.maximum(
        np.log(lim["t"][use] / lim_shares.t_max[use]),
        np.log(lim_shares.t_min[use] / lim["t"][use]),
    )
    return None


def gate(args: argparse.Namespace) -> int:
    """Compare the limited and the unlimited run at equal snapshot times; return the exit class."""
    rep = attenuation.Report()
    params = check_premises(args)
    first_path = os.path.join(args.limited, "snap", "snapshot_0001.hdf5")
    first = attenuation.load_snapshot(first_path, FIELDS)
    t_init = float(params.get("SPH:initial_temperature", "nan"))
    attenuation.require_finite(
        "first snapshot", Temperatures=first["Temperatures"], t_init=np.array([t_init])
    )
    t_scale = t_init / float(np.median(first["Temperatures"]))
    lim_run = Run(args.limited, first, t_scale)
    ctl_run = Run(args.control, first, t_scale)
    cloud = locate_cloud(first, attenuation.star_particles(first_path)[0][0])

    numbers = sorted(
        set(snapshot_numbers(args.limited)) & set(snapshot_numbers(args.control))
    )
    compared: List[int] = []
    saved: Dict[str, np.ndarray] = {}
    ended = None
    for number in [n for n in numbers if n >= 2]:
        lim, ctl = lim_run.get(number), ctl_run.get(number)
        if not np.isclose(lim["time"], ctl["time"], rtol=TIME_MATCH, atol=0.0):
            raise attenuation.GateExit(
                attenuation.EXIT_VOID,
                f"void run, snapshot {number} has the time {lim['time']:.9e} in the limited "
                f"run and {ctl['time']:.9e} in the control run",
            )
        both = cloud.domain & lim["lit"] & ctl["lit"]
        if int(both.sum()) < MIN_COMPARE:
            continue
        neighbours = upstream_neighbours(ctl, cloud.star, cloud.domain, cloud.dist)
        ctl_shares = share_bounds(ctl_run, number, neighbours, cloud.dist)
        ended = control_leaves_range(rep, number, ctl, ctl_shares, cloud)
        if ended is not None:
            break
        lim_shares = share_bounds(lim_run, number, neighbours, cloud.dist)
        ended = compare_snapshot(
            rep,
            number,
            (lim, ctl),
            (lim_shares, ctl_shares),
            both,
            cloud,
            not compared,
            saved,
        )
        if ended is not None:
            break
        compared.append(number)
    if ended is not None:
        rep.info("end of comparison", ended)
    elif compared:
        rep.info(
            "end of comparison",
            f"the last snapshot both runs hold ({numbers[-1]}) is reached; the control run "
            "stays inside the range of its own shares and the H2 stays frozen",
        )
    if compared:
        rep.info(
            "window",
            f"snapshots {compared[0]} to {compared[-1]} ({len(compared)} read); "
            "the snapshots before it have fewer than "
            f"{MIN_COMPARE} particles lit in both runs",
        )
    elif not rep.failed:
        raise attenuation.GateExit(
            attenuation.EXIT_VOID,
            f"void run, no snapshot has {MIN_COMPARE} particles lit in both runs",
        )
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
    except attenuation.GateExit as stop:
        print(f"FAIL  {stop.message}")
        if stop.code == attenuation.EXIT_MISSING:
            print(
                "INFO  run ./run.sh first: it writes the control and limited run directories"
            )
        code = stop.code
    print(f"exit class {code}")
    return code


if __name__ == "__main__":
    sys.exit(main())
