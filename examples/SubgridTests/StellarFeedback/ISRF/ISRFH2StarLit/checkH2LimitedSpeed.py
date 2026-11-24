#!/usr/bin/env python3
"""Check the H2 line absorption of the ISRF transport at a limited light speed.

Both checks read runs of the star-lit box (`ISRFH2StarLit`) with `ISRF_packet_light_speed: cap`,
a non-periodic box and `--temperature` instead of `--cooling`. The cooling is initialised but never
called, so the species and the temperature stay exactly at their initial values. The light then
settles to a steady state once its front has passed, whatever its delay.

Subcommands
-----------
limited   Line transmission of the steady light against the shielding function at the straight-ray
          column (`--run DIR`, the directory of one leg).
ripple    Time structure of the light of a star that steps coarsely (`--coarse DIR --fine DIR`).

Exit codes
----------
0 all checks passed
1 a bar failed
2 a required file, field or key is missing
3 a value is not finite
4 the run is void: a premise or a positive control is not met

The reference of `limited` is the Wolcott-Green and Haiman (2019) fit of `checkH2Attenuation.py`,
evaluated at the column n_H2 r of each particle. It does not use the transport.

Bar of `limited`. In frozen gas the carried column is a sum of the mean n_H2 of the two ends of each
hop times the hop path, and the paths add up to r. Every share that a particle holds, whatever its
delay, then has a column between n_min r and n_max r, with n_min and n_max the extreme H2 densities
of the read domain. The pooled transmission is a mean of the shares' transmissions, so it lies
between the fit at the two columns. The bar of a particle is the larger of the two logarithmic
distances, plus the float rounding of the fit (as in `checkH2Attenuation.py starlit`), plus the
largest change of its own ratio over the last two snapshot intervals (the light is steady only to
that change). The delay of the light does not enter the bar. It enters one premise: every read
particle is lit, and its ratio changes by less than its own column bar over the last two intervals.

Bar of `ripple`. This is a differential test: the reference is the same box with a star that steps
every gas step, run with the same binary. It sees only faults that depend on the star step. A star
that steps every 8 gas steps must give the same mean light as one that steps every gas step, and
no modulation at the star-step period. Three lines, each on shell means of the own LW energy over
the last 16 snapshots (two star steps):
- Energy: the time mean of the coarse run over that of the fine run must be 1. Energy
  conservation makes the long-time ratio 1. The bar is the peak-to-peak range of the fine run
  over its window mean: the fine run shows how far a window mean of this box can lie from the
  long-time mean. The range of the coarse run is not in the bar, so the run under test cannot
  widen its own tolerance.
- Noise: the peak-to-peak range of the coarse run must not exceed that of the fine run. A coarse
  star adds light at the star-step period only, and the ripple line gates that period. Any further
  variation in the coarse run is not explained by the star step.
- Ripple: the amplitude of the Fourier component of the coarse series at the star-step period over
  two periods must stay under the half peak-to-peak range of the fine run in the same shell.
The noise and ripple bars are a floor read from the control run, provisional until the operator
rules on card h2_flimit_coarse_bars_1: no mechanism gives the size of the variation of a correct
run.

Both legs must come from the same binary: the Revision, the configure line and the geajay
revision of the two start-up logs are compared and a difference voids the run. `limited` prints
the revisions it finds and cannot check them against an expected value.

Far-light warning. The start-up log of `cap_coarse` and `cap_half` prints that more than 3 per cent
of the light was booked in far arrival cells. That limit was chosen for light that crosses a
smoothing length in 2 to 3 gas steps. At 0.99 c the light needs about 2.8 gas steps, at 0.5 c about
5.6, and a coarse star books its light over 8 gas steps, so these legs lie outside that range. The
warning does not change the steady transmission that `limited` reads.
"""

import argparse
import math
import os
import re
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import checkH2Attenuation as base  # noqa: E402

EXIT_OK, EXIT_BAR, EXIT_MISSING, EXIT_NONFINITE, EXIT_VOID = 0, 1, 2, 3, 4
GateExit = base.GateExit
Report = base.Report

UNIT_PC_CGS = 3.0856775814913673e18
HYDROGEN_MASS_FRACTION = 0.76
H2_FLOAT_REL = 2.0**-23
PLANT_SIZE = 0.05
STAR_PERIOD_STEPS = 8
STAR_STEP_KEY = "Stars:max_timestep_young_Myr"
WINDOW = 2 * STAR_PERIOD_STEPS
SHELL_WIDTH_IN_H = 0.5
REQUIRED_SERIES = (
    "ParticleIDs",
    "Coordinates",
    "SmoothingLengths",
    "Densities",
    "Masses",
)


def snapshot_path(run: str, index: int) -> str:
    return os.path.join(run, "snap", f"snapshot_{index:04d}.hdf5")


def last_snapshot_index(run: str) -> int:
    """Return the highest snapshot number of a leg directory, void if there is none."""
    folder = os.path.join(run, "snap")
    try:
        names = [n for n in os.listdir(folder) if n.startswith("snapshot_")]
    except OSError as err:
        raise GateExit(EXIT_MISSING, f"no snapshot directory in {run}: {err}")
    if not names:
        raise GateExit(EXIT_MISSING, f"no snapshot in {folder}")
    return max(int(n[len("snapshot_") : -len(".hdf5")]) for n in names)


def premise_parameters(params: Dict[str, str], rep: Report, log_path: str) -> None:
    """Void the run unless the parameters are those of a limited-speed, frozen-species leg."""
    wanted = [
        ("GEARFeedback:ISRF_H2_line_absorption", "1", None),
        ("GEARFeedback:with_photoionization", "0", "0"),
        ("GEARFeedback:HII_couple_ionization_rate", "0", "0"),
        ("GEARFeedback:with_stellar_wind_feedback", "0", None),
        ("GEARFeedback:supernovae_efficiency", "0", None),
        ("GEARFeedback:ISRF_packet_light_speed", "cap", None),
        ("InitialConditions:periodic", "0", None),
    ]
    bad = [k for k, v, d in wanted if not base.param_equals(params, k, v, d)]
    if bad:
        raise GateExit(
            EXIT_VOID, f"void run, premise parameters not as required: {bad}"
        )
    try:
        text = open(log_path).read()
    except OSError as err:
        raise GateExit(EXIT_MISSING, f"start-up log {log_path} is missing: {err}")
    if not re.search(r"h2_transported_lw_rate_shielding = 1\b", text):
        raise GateExit(
            EXIT_VOID,
            "void run, the start-up line h2_transported_lw_rate_shielding = 1 is absent",
        )
    rep.info(
        "premise parameters",
        "HII, winds, SN off, cap light speed, isolated box, line absorption on",
    )


BUILD_LINES = (
    ("SWIFT revision", r"^\s*Revision:.*$"),
    ("configure line", r"^\s*Config\. options:.*$"),
    ("geajay revision", r"geajay version = .*$"),
)


def read_build(log_path: str) -> Dict[str, Optional[str]]:
    """Return the build lines of a start-up log; a line that the log lacks is None."""
    try:
        text = open(log_path).read()
    except OSError as err:
        raise GateExit(EXIT_MISSING, f"start-up log {log_path} is missing: {err}")
    out: Dict[str, Optional[str]] = {}
    for name, pattern in BUILD_LINES:
        found = re.search(pattern, text, re.MULTILINE)
        out[name] = found.group(0).strip() if found else None
    return out


def premise_same_binary(log_a: str, log_b: str, rep: Report) -> None:
    """Void the run when a build line that both logs hold differs; print what cannot be compared."""
    a, b = read_build(log_a), read_build(log_b)
    for name, _ in BUILD_LINES:
        if a[name] is None or b[name] is None:
            rep.info("binary", f"{name} is not in both logs: it cannot be checked")
        elif a[name] != b[name]:
            raise GateExit(
                EXIT_VOID,
                f"void run, the two legs differ in the {name}: {a[name]!r} and {b[name]!r}",
            )
        else:
            rep.info("binary", f"{name} equal in both legs: {a[name][:100]}")


def print_build(log_path: str, rep: Report) -> None:
    """Print the build lines of one log; one leg has no expected value to compare with."""
    for name, line in read_build(log_path).items():
        shown = line[:100] if line else "not in the log"
        rep.info(
            "binary", f"{name}: {shown} (printed only, not checked against a reference)"
        )


def premise_frozen(first: dict, last: dict, h2_ratio: float, rep: Report):
    """Void the run unless H2 and the internal energy are bit-equal and H2 is the start value."""
    if not np.array_equal(first["ParticleIDs"], last["ParticleIDs"]):
        raise GateExit(
            EXIT_VOID, "void run, the two snapshots hold different particle orders"
        )
    if not np.array_equal(first["H2I"], last["H2I"]):
        n = int(np.sum(first["H2I"] != last["H2I"]))
        raise GateExit(
            EXIT_VOID, f"void run, H2I changed on {n} particles: the cooling ran"
        )
    if not np.array_equal(first["InternalEnergies"], last["InternalEnergies"]):
        raise GateExit(
            EXIT_VOID, "void run, the internal energy changed: the cooling ran"
        )
    start = 2.0 * h2_ratio * HYDROGEN_MASS_FRACTION
    off = float(np.max(np.abs(last["H2I"] / start - 1.0)))
    if off > H2_FLOAT_REL:
        raise GateExit(
            EXIT_VOID,
            f"void run, H2I differs from the start value {start:.6g} by {off:.3g} (float rel "
            f"{H2_FLOAT_REL:.3g}): --temperature did not start the species",
        )
    rep.info(
        "frozen species", f"H2I and u bit-equal first to last; H2I equals {start:.6g}"
    )


def ratio_series(path_list: List[str]) -> List[Tuple[np.ndarray, np.ndarray]]:
    """Return (own energy, own dissociating energy) of each snapshot, float64."""
    out = []
    for p in path_list:
        s = base.load_snapshot(p, ("ParticleIDs",) + base.REQUIRED_OWN)
        base.require_finite(p, **{k: s[k] for k in base.REQUIRED_OWN})
        out.append((s["LWOwnSpecificEnergies"], s["LWOwnDissociatingSpecificEnergies"]))
    return out


def own_ratio(e: np.ndarray, d: np.ndarray) -> np.ndarray:
    return d / np.where(e > 0.0, e, 1.0)


def unsteady_term(
    series: List[Tuple[np.ndarray, np.ndarray]], lit: np.ndarray
) -> np.ndarray:
    """Return per particle the largest |ln| change of the own ratio over the given snapshots."""
    ratios = [own_ratio(e, d) for e, d in series]
    if any(np.any(r[lit] <= 0.0) for r in ratios):
        raise GateExit(
            EXIT_VOID,
            "void run, a read particle has no dissociating energy in one of the last three snapshots "
            "(not yet lit, or the light is fully absorbed)",
        )
    out = np.zeros(lit.shape)
    for a, b in zip(ratios[:-1], ratios[1:]):
        out = np.maximum(
            out, np.abs(np.log(np.where(lit, b / np.where(lit, a, 1.0), 1.0)))
        )
    return out


def smallest_detectable(
    dev: np.ndarray, bar: np.ndarray, lit: np.ndarray
) -> Tuple[float, float]:
    """Return the smallest uniform ln error, up and down, that puts one lit probe over its bar."""
    up = float(np.min((bar - dev)[lit]))
    down = float(np.min((bar + dev)[lit]))
    return up, down


def gate_limited(args: argparse.Namespace) -> int:
    """Check the steady line transmission of one limited-speed leg.

    Premises (each a void run): the parameters of the leg; H2 and the internal energy bit-equal
    between the first and the last snapshot and H2 at the start value of the box; every particle of
    the read domain lit at the last snapshot; the own ratio changes by less than the column bar of
    the particle over the last two intervals; a uniform error of 5 per cent is seen by at least one
    probe (the positive control of the bar).

    Bar. See the module docstring. A check line fails when one lit probe lies outside its bar.
    """
    rep = Report()
    run = args.run
    last_index = args.last if args.last is not None else last_snapshot_index(run)
    if last_index < 3:
        raise GateExit(
            EXIT_VOID, f"void run, the last snapshot {last_index} is too early"
        )
    fix = base.load_json(os.path.join(run, "fixture.json"), "geometry file")
    first = base.load_snapshot(snapshot_path(run, 1), base.REQUIRED_GAS)
    fields = base.REQUIRED_GAS + base.REQUIRED_OWN
    last = base.load_snapshot(snapshot_path(run, last_index), fields)
    base.require_finite("last snapshot", **{k: last[k] for k in fields})
    base.require_finite("first snapshot", **{k: first[k] for k in base.REQUIRED_GAS})
    base.check_box_faces(first)
    params = base.read_params(os.path.join(run, "used_parameters.yml"))
    premise_parameters(params, rep, os.path.join(run, "output.log"))
    print_build(os.path.join(run, "output.log"), rep)
    base.calibrate_temperatures(first, last, params, rep)

    r, h, mask = base.probe_geometry(last, np.array(fix["star_kpc"]))
    premise_frozen(first, last, float(fix["h2_ratio"]), rep)
    n_h2 = base.h2_column_density(last)
    own, meas = base.own_transmission(last)
    lit = mask & (own > 0.0)
    if int(lit.sum()) != int(mask.sum()):
        raise GateExit(
            EXIT_VOID,
            f"void run, {int(mask.sum() - lit.sum())} of {int(mask.sum())} read particles are "
            "not lit: the light front has not passed",
        )
    series = ratio_series(
        [snapshot_path(run, k) for k in range(last_index - 2, last_index + 1)]
    )
    unsteady = unsteady_term(series, lit)
    if args.plant_t_scale is not None:
        e32, d32 = series[-1]
        d_plant = (d32.astype(np.float32) * np.float32(args.plant_t_scale)).astype(
            np.float64
        )
        meas = np.where(own > 0.0, d_plant / np.where(own > 0.0, own, 1.0), meas)
        rep.info("plant", f"own dissociating energy times {args.plant_t_scale:g}")
    base.require_finite("own ratio", ratio=meas[lit])

    dens = base.gas_number_density(last)
    temp = last["Temperatures"]
    column = n_h2 * r
    exp = base.transmission(column, temp, dens)
    slope = base.log_slope(np.maximum(column, 1.0), temp, dens)
    n_both = np.concatenate([n_h2[mask], base.h2_column_density(first)[mask]])
    base.require_finite("n_H2 of the read domain", n_h2=n_both)
    n_range = (float(np.min(n_both)), float(np.max(n_both)))
    float_term = (
        (r / np.maximum(h, 1e-300) + 2.0) * (2.0 + 2.0 * slope) * base.FLOAT32_EPS
    )
    column_bar = base.column_spread_bar(r, n_h2, mask, temp, dens, n_range) + float_term
    base.require_finite("column bar", bar=column_bar[lit], unsteady=unsteady[lit])
    too_unsteady = lit & (unsteady > column_bar)
    if too_unsteady.any():
        raise GateExit(
            EXIT_VOID,
            f"void run, the own ratio of {int(too_unsteady.sum())} particles changed by more than "
            f"their column bar over the last two intervals (largest {float(unsteady.max()):.3g})",
        )
    bar = column_bar + unsteady
    dev = np.log(np.maximum(meas, 1e-300) / exp)
    base.require_finite("deviation", dev=dev[lit], bar=bar[lit])
    over = lit & (np.abs(dev) > bar)
    up, down = smallest_detectable(dev, bar, lit)
    rep.info(
        "column spread",
        f"n_H2 over the domain {n_range[0]:.4g} to {n_range[1]:.4g} cm^-3; "
        f"median bar {float(np.median(bar[lit])):.3g}, median unsteady term "
        f"{float(np.median(unsteady[lit])):.3g}",
    )
    rep.info(
        "smallest detectable uniform error",
        f"ln factor +{up:.4g} (transmission too high) or -{down:.4g} (too low), "
        f"i.e. a factor {math.exp(up):.4f} or {math.exp(-down):.4f}",
    )
    sees_both = up < math.log(1.0 + PLANT_SIZE) and down < -math.log(1.0 - PLANT_SIZE)
    if not sees_both and args.plant_t_scale is None:
        raise GateExit(
            EXIT_VOID,
            "void run, the bar cannot see a uniform error of 5 per cent in both signs",
        )
    rep.check(
        "own transmission against f_sh(n_H2 r)",
        not over.any(),
        f"{int(over.sum())} of {int(lit.sum())} probes outside |ln(own ratio/f)| <= bar; "
        f"worst margin {float(np.max((np.abs(dev) / bar)[lit])):.3g} of the bar; "
        f"median ln(own ratio/f) {float(np.median(dev[lit])):.4g}",
    )
    return rep.exit_code()


def load_series(
    run: str, indices: List[int]
) -> Tuple[Dict[str, np.ndarray], np.ndarray, np.ndarray]:
    """Return the static fields, the own energies (snapshots x particles) and the snapshot times."""
    first = base.load_snapshot(
        snapshot_path(run, indices[0]), REQUIRED_SERIES + ("LWOwnSpecificEnergies",)
    )
    base.require_finite(
        f"{run} static fields", **{k: first[k] for k in REQUIRED_SERIES[1:]}
    )
    energies, times = [], []
    for k in indices:
        s = base.load_snapshot(
            snapshot_path(run, k), ("ParticleIDs", "LWOwnSpecificEnergies")
        )
        if not np.array_equal(s["ParticleIDs"], first["ParticleIDs"]):
            raise GateExit(
                EXIT_VOID, f"void run, snapshot {k} of {run} has another particle order"
            )
        base.require_finite(f"{run} snapshot {k}", own=s["LWOwnSpecificEnergies"])
        energies.append(s["LWOwnSpecificEnergies"])
        times.append(s["time"])
    return first, np.array(energies), np.array(times)


def read_updates(path: str) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return the step times and the hydro (gas) and star update counts of timesteps.txt."""
    rows = []
    try:
        with open(path) as f:
            for line in f:
                tok = line.split()
                if tok and not tok[0].startswith("#") and len(tok) >= 10:
                    rows.append((float(tok[1]), int(tok[7]), int(tok[9])))
    except (OSError, ValueError) as err:
        raise GateExit(EXIT_MISSING, f"timesteps file {path} is unreadable: {err}")
    if not rows:
        raise GateExit(EXIT_MISSING, f"timesteps file {path} holds no step")
    a = np.array(rows)
    return a[:, 0], a[:, 1].astype(int), a[:, 2].astype(int)


def premise_same_runs(
    coarse: Dict[str, str], fine: Dict[str, str], rep: Report
) -> None:
    """Void the run unless the parameters of the two legs differ only in the star step."""
    differing = sorted(
        k for k in set(coarse) | set(fine) if coarse.get(k) != fine.get(k)
    )
    if differing != [STAR_STEP_KEY]:
        raise GateExit(
            EXIT_VOID,
            f"void run, the two legs differ in {differing}, not only in {STAR_STEP_KEY}",
        )
    rep.info("parameters", f"the legs differ only in {STAR_STEP_KEY}")


def premise_steps(
    run: str, n_gas: int, star_every: int, indices: List[int], rep: Report
) -> None:
    """Void the run unless every gas particle steps at each row and the star steps as stated."""
    _, gas, star = read_updates(os.path.join(run, "timesteps.txt"))
    rows = slice(indices[0], indices[-1] + 1)
    if not np.all(gas[rows] == n_gas):
        raise GateExit(
            EXIT_VOID, f"void run, not every gas particle steps in each step of {run}"
        )
    spaces = np.diff(np.flatnonzero(star[rows] > 0))
    ok = len(spaces) >= 1 and np.all(spaces == star_every)
    if not ok:
        raise GateExit(
            EXIT_VOID,
            f"void run, the star of {run} does not step every {star_every} gas steps in the window "
            f"(gaps {sorted(set(spaces.tolist()))})",
        )
    rep.info(
        "steps", f"{run}: all gas steps together, star steps every {star_every} rows"
    )


def shell_series(
    static: Dict[str, np.ndarray], energies: np.ndarray, star: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Return the volume-weighted shell means of the own energy and the shell edges in pc."""
    box, ul = static["box"], static["unit_length"]
    d = static["Coordinates"] - star[None, :]
    r = np.linalg.norm(d, axis=1) * ul / UNIT_PC_CGS
    h = static["SmoothingLengths"]
    h_typical = float(np.median(h))
    r_min = 3.0 * base.KERNEL_GAMMA * h_typical * ul / UNIT_PC_CGS
    r_max = (0.5 * box - 2.0 * h_typical) * ul / UNIT_PC_CGS
    width = SHELL_WIDTH_IN_H * h_typical * ul / UNIT_PC_CGS
    edges = np.arange(r_min, r_max, width)
    volume = static["Masses"] / static["Densities"]
    out = np.zeros((energies.shape[0], len(edges) - 1))
    for j, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        sel = (r >= lo) & (r < hi)
        if not sel.any():
            raise GateExit(
                EXIT_VOID, f"void run, the shell {lo:.3g} to {hi:.3g} pc holds no gas"
            )
        if np.any(energies[:, sel] <= 0.0):
            raise GateExit(
                EXIT_VOID,
                f"void run, a particle of the shell {lo:.3g} to {hi:.3g} pc is not lit in "
                "every snapshot of the window",
            )
        out[:, j] = (energies[:, sel] * volume[sel]).sum(axis=1) / volume[sel].sum()
    return out, edges


def peak_to_peak(series: np.ndarray) -> np.ndarray:
    return (series.max(axis=0) - series.min(axis=0)) / series.mean(axis=0)


def period_amplitude(series: np.ndarray, period: int) -> np.ndarray:
    """Return the Fourier amplitude per shell at a period that divides the window."""
    n = series.shape[0]
    if n % period:
        raise GateExit(
            EXIT_VOID, f"{n} snapshots are not a whole number of {period}-step periods"
        )
    phase = np.exp(-2j * np.pi * np.arange(n) / period)
    return (
        np.abs(((series / series.mean(axis=0)) * phase[:, None]).sum(axis=0)) * 2.0 / n
    )


def plant_ripple(
    energies: np.ndarray,
    indices: List[int],
    r_pc: np.ndarray,
    amplitude: float,
    step_pc: float,
) -> np.ndarray:
    """Return the energies times 1 + A cos(2 pi (k - r / step) / period), a wave of the star step."""
    k = np.asarray(indices, float)[:, None]
    phase = 2.0 * np.pi * (k - r_pc[None, :] / step_pc) / STAR_PERIOD_STEPS
    return (
        (energies * (1.0 + amplitude * np.cos(phase)))
        .astype(np.float32)
        .astype(np.float64)
    )


def gate_ripple(args: argparse.Namespace) -> int:
    """Check the light of a coarse star against the same box with a star that steps every gas step.

    Premises (each a void run): the parameters of the two legs differ in nothing but the star
    step; the logs show the same SWIFT revision, configure line and geajay revision; H2 and the
    internal energy of each leg are bit-equal between its first and last snapshot; the snapshots
    hold the same particles at the same times; every gas particle steps in every step; the star of
    the coarse leg steps every 8 gas steps and that of the fine leg every step; every particle of
    the shells is lit in every snapshot of the window; every value read is finite.
    """
    rep = Report()
    coarse_dir, fine_dir = args.coarse, args.fine
    n_last = last_snapshot_index(coarse_dir)
    if n_last != last_snapshot_index(fine_dir):
        raise GateExit(EXIT_VOID, "void run, the legs end at different snapshots")
    indices = list(range(n_last - WINDOW + 1, n_last + 1))
    fix = base.load_json(os.path.join(coarse_dir, "fixture.json"), "geometry file")
    all_params = {}
    for d in (coarse_dir, fine_dir):
        all_params[d] = base.read_params(os.path.join(d, "used_parameters.yml"))
        premise_parameters(all_params[d], rep, os.path.join(d, "output.log"))
    premise_same_runs(all_params[coarse_dir], all_params[fine_dir], rep)
    premise_same_binary(
        os.path.join(coarse_dir, "output.log"),
        os.path.join(fine_dir, "output.log"),
        rep,
    )
    for d in (coarse_dir, fine_dir):
        first = base.load_snapshot(snapshot_path(d, 1), base.REQUIRED_GAS)
        last = base.load_snapshot(snapshot_path(d, n_last), base.REQUIRED_GAS)
        base.require_finite(
            f"{d} first snapshot", **{k: first[k] for k in base.REQUIRED_GAS}
        )
        base.require_finite(
            f"{d} last snapshot", **{k: last[k] for k in base.REQUIRED_GAS}
        )
        premise_frozen(first, last, float(fix["h2_ratio"]), rep)
    static_c, e_c, t_c = load_series(coarse_dir, indices)
    static_f, e_f, t_f = load_series(fine_dir, indices)
    if not np.allclose(t_c, t_f, rtol=base.TIME_RTOL, atol=0.0):
        raise GateExit(EXIT_VOID, "void run, the snapshot times of the two legs differ")
    n_gas = len(static_c["ParticleIDs"])
    premise_steps(coarse_dir, n_gas, STAR_PERIOD_STEPS, indices, rep)
    premise_steps(fine_dir, n_gas, 1, indices, rep)
    star = np.array(fix["star_kpc"])
    step_pc = float(args.step_pc)
    if args.plant_coarse_energy is not None:
        e_c = (
            (e_c * np.float32(args.plant_coarse_energy))
            .astype(np.float32)
            .astype(np.float64)
        )
        rep.info("plant", f"coarse own energy times {args.plant_coarse_energy:g}")
    if args.plant_ripple is not None:
        r_pc = np.linalg.norm(static_c["Coordinates"] - star[None, :], axis=1)
        r_pc = r_pc * static_c["unit_length"] / UNIT_PC_CGS
        e_c = plant_ripple(e_c, indices, r_pc, args.plant_ripple, step_pc)
        rep.info(
            "plant",
            f"coarse own energy modulated by {args.plant_ripple:g} at the star period",
        )
    s_c, edges = shell_series(static_c, e_c, star)
    s_f, _ = shell_series(static_f, e_f, star)
    base.require_finite("shell series", coarse=s_c, fine=s_f)
    if np.any(s_c <= 0.0) or np.any(s_f <= 0.0):
        raise GateExit(EXIT_VOID, "void run, a shell holds no light in the window")
    rep.info(
        "shells",
        f"{len(edges) - 1} shells from {edges[0]:.3g} to {edges[-1]:.3g} pc, "
        f"window {WINDOW} snapshots",
    )

    ratio = np.log(s_c.mean(axis=0) / s_f.mean(axis=0))
    bar_energy = peak_to_peak(s_f)
    noise_coarse, noise_fine = peak_to_peak(s_c), peak_to_peak(s_f)
    amp = period_amplitude(s_c, STAR_PERIOD_STEPS)
    bar_ripple = 0.5 * peak_to_peak(s_f)
    base.require_finite(
        "bars", energy=bar_energy, ripple=bar_ripple, noise=noise_coarse, amp=amp
    )
    rep.info(
        "smallest detectable",
        f"energy error {float(bar_energy.min()):.3g} in ln "
        f"(largest shell {float(bar_energy.max()):.3g}); "
        f"ripple amplitude {float(bar_ripple.min()):.3g} "
        f"(largest shell {float(bar_ripple.max()):.3g}); "
        "the ripple and noise bars are a floor from the control run, provisional",
    )
    rep.check(
        "time-mean energy, coarse over fine",
        bool(np.all(np.abs(ratio) <= bar_energy)),
        f"largest |ln ratio| {float(np.max(np.abs(ratio))):.4g}, margin "
        f"{float(np.max(np.abs(ratio) / bar_energy)):.3g} of the bar; ln ratio by shell "
        f"{np.round(ratio, 4).tolist()}",
    )
    rep.check(
        "coarse leg no noisier than the fine leg (floor from the control run, provisional)",
        bool(np.all(noise_coarse <= noise_fine)),
        f"largest range of the coarse leg {float(noise_coarse.max()):.3g}, of the fine leg "
        f"{float(noise_fine.max()):.3g}; worst ratio "
        f"{float(np.max(noise_coarse / noise_fine)):.3g}",
    )
    rep.check(
        "no ripple at the star-step period (floor from the control run, provisional)",
        bool(np.all(amp <= bar_ripple)),
        f"largest amplitude {float(amp.max()):.4g}, "
        f"margin {float(np.max(amp / bar_ripple)):.3g} of the bar; fine-run amplitude at the same "
        f"period {float(period_amplitude(s_f, STAR_PERIOD_STEPS).max()):.3g}",
    )
    return rep.exit_code()


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter
    )
    sub = p.add_subparsers(dest="cmd", required=True)
    sp = sub.add_parser("limited", help="steady line transmission of one leg")
    sp.add_argument(
        "--run", required=True, help="leg directory (snap/, used_parameters.yml, ...)"
    )
    sp.add_argument(
        "--last",
        type=int,
        default=None,
        help="last snapshot number, default the highest",
    )
    sp.add_argument(
        "--plant-t-scale",
        type=float,
        default=None,
        help="offline plant: own dissociating energy times this",
    )
    sp.set_defaults(func=gate_limited)
    sp = sub.add_parser("ripple", help="coarse-star light against the fine-star light")
    sp.add_argument("--coarse", required=True)
    sp.add_argument("--fine", required=True)
    sp.add_argument(
        "--step-pc",
        type=float,
        default=0.0696,
        help="light travel per gas step in pc (phase of a planted ripple)",
    )
    sp.add_argument(
        "--plant-coarse-energy",
        type=float,
        default=None,
        help="offline plant: coarse own energy times this",
    )
    sp.add_argument(
        "--plant-ripple",
        type=float,
        default=None,
        help="offline plant: modulation amplitude at the star period",
    )
    sp.set_defaults(func=gate_ripple)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except GateExit as exc:
        print(f"{'VOID' if exc.code == EXIT_VOID else 'STOP'}  {exc.message}")
        return exc.code


if __name__ == "__main__":
    sys.exit(main())
