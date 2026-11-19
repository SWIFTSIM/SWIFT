#!/usr/bin/env python3
"""Check the H2 line absorption of the ISRF transport on the three H2 fixtures.

Subcommands
-----------
starlit     Transmission at depth against the shielding function (star-lit box).
twosource   Merge of two lights of different columns (two sources).
velocity    Reset of the carried column by a velocity shift.
parity      Per-particle comparison of two runs (two ranks, or restart).
refusal     A run that must be refused at start-up was refused.

Exit codes
----------
0 all checks passed
1 a bar failed
2 a required file, field or key is missing
3 a value is not finite
4 the run is void: a premise or a positive control is not met

A missing field, a void run and a non-finite value never pass and never crash.

The shielding function is the fit of the design: x = N_H2 / 5e14 cm^-2,
b_5 = sqrt(2 k T / m_H2) / (1 km/s), f = 0.965 / (1 + x / b_5)^alpha +
0.035 exp(-8.5e-4 sqrt(1 + x)) / sqrt(1 + x), with alpha from the temperature
and density fit of Wolcott-Green and Haiman (2019). The density is floored
at 1e-8 cm^-3 as in the module. The expected transmission is f(N) / f(0).
"""

import argparse
import json
import math
import re
import sys
from typing import Dict, List, Optional, Tuple

import h5py
import numpy as np

EXIT_OK, EXIT_BAR, EXIT_MISSING, EXIT_NONFINITE, EXIT_VOID = 0, 1, 2, 3, 4

K_B = 1.380649e-16
M_H = 1.6735575e-24
FLOAT32_EPS = 2.0**-24
MERGE_POINT_RTOL = 1e-5
FROZEN_H2_REL = 1e-3
TEMPERATURE_REL = 0.05
THIN_COLUMN_CM2 = 1e17
WIDTH_RESET = 3.0
KERNEL_GAMMA = 1.9364916731  # Wendland C2 support radius in units of h
PARITY_COLUMN = "LWLineH2ColumnDensities"

SPECIES_MASS_NUMBER = {
    "HI": 1.0,
    "HII": 1.0,
    "HM": 1.0,
    "HeI": 4.0,
    "HeII": 4.0,
    "HeIII": 4.0,
    "H2I": 2.0,
    "H2II": 2.0,
    "DI": 2.0,
    "DII": 2.0,
    "HDI": 3.0,
}
ELECTRONS_PER_MASS = {
    "HII": 1.0,
    "HeII": 0.25,
    "HeIII": 0.5,
    "H2II": 0.5,
    "DII": 0.5,
    "HM": -1.0,
}
REQUIRED_LW = (
    "LWSpecificEnergies",
    "LWDissociatingSpecificEnergies",
    "LWLineTransmissions",
)
REQUIRED_OWN = (
    "LWOwnSpecificEnergies",
    "LWOwnDissociatingSpecificEnergies",
)
REQUIRED_GAS = (
    "ParticleIDs",
    "Coordinates",
    "SmoothingLengths",
    "Densities",
    "InternalEnergies",
    "H2I",
) + tuple(SPECIES_MASS_NUMBER)
FACE_MARGIN = 1e-6


class GateExit(Exception):
    """Stop the gate with an exit class and a message."""

    def __init__(self, code: int, message: str):
        super().__init__(message)
        self.code = code
        self.message = message


class Report:
    """Collect named checks and print one line for each."""

    def __init__(self) -> None:
        self.failed: List[str] = []

    def check(self, name: str, ok: bool, detail: str) -> None:
        print(f"{'PASS' if ok else 'FAIL'}  {name}: {detail}")
        if not ok:
            self.failed.append(name)

    @staticmethod
    def info(name: str, detail: str) -> None:
        print(f"INFO  {name}: {detail}")

    @staticmethod
    def reading(name: str, detail: str) -> None:
        print(f"READ  {name}: {detail}  (printed only, no exit code)")

    def exit_code(self) -> int:
        return EXIT_BAR if self.failed else EXIT_OK


def shielding(n_h2: np.ndarray, temp: np.ndarray, dens_cm3: np.ndarray) -> np.ndarray:
    """Return the H2 line shielding function f_sh(N) of the design.

    Parameters
    ----------
    n_h2 : np.ndarray
        H2 column in cm^-2.
    temp : np.ndarray
        Gas temperature in K.
    dens_cm3 : np.ndarray
        Gas density in cm^-3 (H atoms), floored at 1e-8.
    """
    n_h2, temp, dens_cm3 = np.broadcast_arrays(
        np.asarray(n_h2, float), np.asarray(temp, float), np.asarray(dens_cm3, float)
    )
    x = n_h2 / 5e14
    b5 = np.sqrt(2.0 * K_B * temp / (2.0 * M_H)) / 1e5
    tp = np.log10(np.clip(temp, 100.0, 8000.0))
    npr = np.log10(np.clip(dens_cm3, 1e-8, 1e7))
    alpha = (0.8711 * tp - 1.928) * np.exp(-0.2856 * npr) + (-0.9639 * tp + 3.892)
    f = 0.965 / (1.0 + x / b5) ** alpha + 0.035 * np.exp(
        -8.5e-4 * np.sqrt(1.0 + x)
    ) / np.sqrt(1.0 + x)
    return np.minimum(f, 1.0)


def transmission(n_h2, temp, dens_cm3):
    """Return f_sh(N) / f_sh(0), the expected line transmission."""
    return shielding(n_h2, temp, dens_cm3) / shielding(0.0, temp, dens_cm3)


def log_slope(n_h2, temp, dens_cm3, rel: float = 1e-3):
    """Return |d ln f_sh / d ln N| by a central difference in ln N."""
    n_h2 = np.asarray(n_h2, float)
    hi = shielding(n_h2 * (1.0 + rel), temp, dens_cm3)
    lo = shielding(n_h2 * (1.0 - rel), temp, dens_cm3)
    return np.abs(np.log(hi / lo) / (np.log1p(rel) - np.log1p(-rel)))


def column_for_transmission(target: float, temp: float, dens_cm3: float) -> float:
    """Return the column where the transmission equals target (bisection in ln N)."""
    lo, hi = math.log(1e8), math.log(1e24)
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        if float(transmission(math.exp(mid), temp, dens_cm3)) > target:
            lo = mid
        else:
            hi = mid
    return math.exp(0.5 * (lo + hi))


def load_json(path: Optional[str], what: str) -> dict:
    if path is None:
        raise GateExit(EXIT_MISSING, f"{what} is required and was not given")
    try:
        with open(path) as f:
            return json.load(f)
    except (OSError, ValueError) as err:
        raise GateExit(EXIT_MISSING, f"{what} {path} is missing or unreadable: {err}")


def load_snapshot(path: str, fields: Tuple[str, ...]) -> Dict[str, np.ndarray]:
    """Read the gas fields, the box size and the units of a snapshot.

    Raises
    ------
    GateExit
        With EXIT_MISSING on an absent file or field.
    """
    try:
        f = h5py.File(path, "r")
    except OSError as err:
        raise GateExit(EXIT_MISSING, f"snapshot {path} cannot be opened: {err}")
    with f:
        if "PartType0" not in f:
            raise GateExit(EXIT_MISSING, f"{path} has no PartType0")
        missing = [k for k in fields if k not in f["PartType0"]]
        if missing:
            raise GateExit(EXIT_MISSING, f"{path} lacks the fields {missing}")
        out = {k: np.asarray(f["PartType0"][k][:], dtype=np.float64) for k in fields}
        out["box"] = np.asarray(f["Header"].attrs["BoxSize"], float)[0]
        out["time"] = (
            float(np.asarray(f["Header"].attrs["Time"]).ravel()[0])
            if "Time" in f["Header"].attrs
            else float("nan")
        )
        out["unit_length"] = float(
            np.asarray(f["Units"].attrs["Unit length in cgs (U_L)"]).ravel()[0]
        )
        out["unit_mass"] = float(
            np.asarray(f["Units"].attrs["Unit mass in cgs (U_M)"]).ravel()[0]
        )
        out["unit_time"] = float(
            np.asarray(f["Units"].attrs["Unit time in cgs (U_t)"]).ravel()[0]
        )
    if "InternalEnergies" in out:
        out["Temperatures"] = temperature_from_energy(out)
    return out


def temperature_from_energy(snap: dict) -> np.ndarray:
    """Return T = (2/3) u mu m_H / k_B in K, with mu from the species mass fractions.

    The species fractions of the snapshot sum to slightly more than 1 (H2 is added on top
    of HI), so mu is a first estimate; calibrate_temperatures() scales it to the IC.
    """
    n_per_mass = sum(snap[k] / a for k, a in SPECIES_MASS_NUMBER.items())
    n_per_mass = n_per_mass + sum(c * snap[k] for k, c in ELECTRONS_PER_MASS.items())
    mu = 1.0 / n_per_mass
    u_cgs = snap["InternalEnergies"] * (snap["unit_length"] / snap["unit_time"]) ** 2
    return (2.0 / 3.0) * u_cgs * mu * M_H / K_B


def calibrate_temperatures(first: dict, last: dict, params: dict, rep: Report) -> None:
    """Scale both temperature arrays so that the first snapshot has SPH:initial_temperature.

    The code sets the initial energy with its own mean molecular weight; the estimate of
    temperature_from_energy() differs from it by a constant factor. The factor cancels in
    the temperature premise, and it is applied to the fit.
    """
    key = "SPH:initial_temperature"
    if key not in params:
        rep.info("temperature scale", f"{key} not in the parameters: factor 1")
        return
    t_init = float(params[key])
    scale = t_init / float(np.median(first["Temperatures"]))
    first["Temperatures"] = first["Temperatures"] * scale
    last["Temperatures"] = last["Temperatures"] * scale
    rep.info(
        "temperature scale",
        f"first snapshot median T estimate scaled by {scale:.5f} to the initial temperature {t_init:g} K",
    )


def check_box_faces(snap: dict, extra: Optional[np.ndarray] = None) -> None:
    """Refuse a snapshot with a coordinate within 1e-6 L of a box face.

    A gas particle on a face can leave the box by less than one ulp of L and be placed in
    a non-existing cell by SWIFT.

    Raises
    ------
    GateExit
        EXIT_VOID with 'coordinate on a box face'.
    """
    box = snap["box"]
    pos = (
        snap["Coordinates"]
        if extra is None
        else np.vstack([snap["Coordinates"], extra])
    )
    d = np.minimum(pos, box - pos)
    bad = int(np.sum(d < FACE_MARGIN * box))
    if bad:
        raise GateExit(
            EXIT_VOID,
            f"coordinate on a box face: {bad} coordinates within {FACE_MARGIN:g} L of a face (fixture IC defect)",
        )


def require_finite(label: str, **arrays: np.ndarray) -> None:
    for name, a in arrays.items():
        if not np.all(np.isfinite(a)):
            raise GateExit(
                EXIT_NONFINITE,
                f"{label}: {name} has {np.sum(~np.isfinite(a))} non-finite values",
            )


def read_params(path: str) -> dict:
    """Return used_parameters.yml as a flat dict of 'Section:key' to value."""
    flat: Dict[str, str] = {}
    section = ""
    try:
        text = open(path).read().splitlines()
    except OSError as err:
        raise GateExit(EXIT_MISSING, f"parameter file {path} is missing: {err}")
    for line in text:
        line = line.split("#")[0].rstrip()
        if not line.strip():
            continue
        if not line.startswith((" ", "\t")):
            section = line.rstrip(":").strip()
            continue
        key, _, val = line.strip().partition(":")
        flat[f"{section}:{key.strip()}"] = val.strip().strip("'\"")
    return flat


def param_equals(params: dict, key: str, expected: str, default: Optional[str]) -> bool:
    val = params.get(key, default)
    if val is None:
        raise GateExit(EXIT_MISSING, f"parameter {key} is absent and has no default")
    try:
        return float(val) == float(expected)
    except ValueError:
        return val == expected


def check_premise_parameters(
    params: dict, rep: Report, mode: str, startup_log: Optional[str]
) -> None:
    """Fail the run as void if the switch-on fixture parameters are not set."""
    wanted = [
        ("GEARFeedback:ISRF_H2_line_absorption", "1", None),
        ("GEARFeedback:with_photoionization", "0", "0"),
        ("GEARFeedback:HII_couple_ionization_rate", "0", "0"),
        ("GEARFeedback:with_stellar_wind_feedback", "0", None),
        ("GEARFeedback:supernovae_efficiency", "0", None),
        ("GEARFeedback:ISRF_packet_light_speed", "unlimited", None),
    ]
    bad = [k for k, v, d in wanted if not param_equals(params, k, v, d)]
    if bad:
        raise GateExit(
            EXIT_VOID, f"void run, premise parameters not as required: {bad}"
        )
    if mode == "a":
        if not param_equals(
            params,
            "GeajayCooling:h2_transported_lw_rate_shielding",
            "background_only",
            None,
        ):
            raise GateExit(
                EXIT_VOID,
                "void run, option a needs h2_transported_lw_rate_shielding: background_only",
            )
    else:
        if startup_log is None:
            raise GateExit(EXIT_MISSING, "option c needs --startup-log")
        text = open(startup_log).read() if startup_log else ""
        if not re.search(r"h2_transported_lw_rate_shielding = 1\b", text):
            raise GateExit(
                EXIT_VOID, "void run, start-up log does not report background_only"
            )
    rep.info(
        "premise parameters",
        "HII, winds, SN off, unlimited speed, library told background_only",
    )


def probe_geometry(
    snap: dict, star_kpc: np.ndarray, h_cm_floor: float = 0.0
) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Return r (cm), h (cm) and the probe mask r < L/2 - 2h around a star."""
    box = snap["box"]
    d = snap["Coordinates"] - star_kpc[None, :]
    d -= box * np.round(d / box)
    cm = snap["unit_length"]
    r = np.linalg.norm(d, axis=1) * cm
    h = snap["SmoothingLengths"] * cm
    mask = (r > 0.0) & (r < 0.5 * box * cm - 2.0 * h)
    return r, h, mask


def gas_number_density(snap: dict) -> np.ndarray:
    """Return rho in H atoms per cm^3 (rho / m_H) in cgs."""
    return snap["Densities"] * snap["unit_mass"] / snap["unit_length"] ** 3 / M_H


def h2_column_density(snap: dict) -> np.ndarray:
    """Return n_H2 = rho X_H2 / (2 m_H) in cm^-3."""
    return (
        snap["Densities"]
        * snap["unit_mass"]
        / snap["unit_length"] ** 3
        * snap["H2I"]
        / (2.0 * M_H)
    )


def probe_premises(
    first: dict,
    last: dict,
    mask: np.ndarray,
    column: np.ndarray,
    rep: Report,
    frozen_h2: bool = True,
) -> np.ndarray:
    """Check frozen H2, the LW source and the temperature; return the read mask.

    With frozen_h2 False the per-probe H2 change is only printed: the caller states its
    own premise on the n_H2 range.

    Raises
    ------
    GateExit
        EXIT_VOID if H2 is not frozen on a read probe or no LW light was emitted.
    """
    if not np.any(last["LWSpecificEnergies"] > 0.0):
        raise GateExit(
            EXIT_VOID,
            "void run, no gas particle holds LW energy (does Z = 0 give LW light?)",
        )
    if not np.array_equal(first["ParticleIDs"], last["ParticleIDs"]):
        order0 = np.argsort(first["ParticleIDs"])
        order1 = np.argsort(last["ParticleIDs"])
        if not np.array_equal(
            first["ParticleIDs"][order0], last["ParticleIDs"][order1]
        ):
            raise GateExit(
                EXIT_VOID, "void run, the two snapshots hold different particles"
            )
        for k in ("H2I", "Temperatures", "Densities"):
            tmp = np.empty_like(first[k])
            tmp[order1] = first[k][order0]
            first[k] = tmp
    dx = np.abs(last["H2I"] - first["H2I"]) / np.maximum(first["H2I"], 1e-300)
    bad_h2 = int(np.sum(mask & (dx > FROZEN_H2_REL)))
    if not frozen_h2:
        rep.info(
            "H2 change per probe",
            f"{bad_h2} probes changed by more than {FROZEN_H2_REL} against the reference snapshot",
        )
    elif bad_h2:
        raise GateExit(
            EXIT_VOID,
            f"void run, H2 changed by more than {FROZEN_H2_REL} on {bad_h2} probes",
        )
    dt = np.abs(last["Temperatures"] - first["Temperatures"]) / first["Temperatures"]
    hot = dt > TEMPERATURE_REL
    skipped = int(np.sum(mask & hot & (column < THIN_COLUMN_CM2)))
    rep.info(
        "temperature premise",
        f"{skipped} thin probes (N < 1e17 cm^-2) skipped, T moved by more than 5 per cent",
    )
    return mask & ~(hot & (column < THIN_COLUMN_CM2))


def hop_length(counters: dict, h_cm: np.ndarray, rep: Report) -> np.ndarray:
    """Return l_hop per particle: the measured mean path, else the smoothing length."""
    if "mean_hop_path_cm" in counters:
        val = float(counters["mean_hop_path_cm"])
        rep.info(
            "l_hop", f"measured mean (d_r - d_s) of the booked shares: {val:.4g} cm"
        )
        return np.full_like(h_cm, val)
    rep.info(
        "l_hop", "packet statistics not given: the smoothing length h of each probe"
    )
    return h_cm


def light_hops(first: dict, last: dict) -> int:
    """Return the number of light hops between the two snapshots, Time_last / Time_first.

    The first snapshot is written after step 1, so the time ratio is the step count of the
    last snapshot when the step is constant (dt_max below the step the particles ask for).

    Raises
    ------
    GateExit
        EXIT_VOID when the ratio is not within 1e-3 of an integer of at least 2.
    """
    ratio = last["time"] / first["time"]
    require_finite("snapshot times", time_ratio=np.array([ratio]))
    n = round(ratio)
    if n < 2 or abs(ratio - n) > 1e-3 * n:
        raise GateExit(
            EXIT_VOID,
            f"void run, time ratio of the snapshots {ratio:.6g} is not a step count of at least 2",
        )
    return int(n)


def own_transmission(last: dict) -> Tuple[np.ndarray, np.ndarray]:
    """Return the own LW field of each particle and its line transmission.

    The own field is the light that crosses the particle itself, before the kernel
    mean. The kernel mean of LWLineTransmissions is a u-weighted mean over the
    neighbours and is biased up near the light front, so the bars read the own ratio.
    """
    own = last["LWOwnSpecificEnergies"]
    ratio = last["LWOwnDissociatingSpecificEnergies"] / np.where(own > 0, own, 1.0)
    return own, ratio


def h2_range_premise(
    first: dict, last: dict, mask: np.ndarray, rep: Report
) -> Tuple[float, float]:
    """Return the n_H2 range of the probe domain over the reference and the last snapshot.

    The carried column is a sum of n_H2 of the receivers at the hop times, so the bar
    needs the range of n_H2 of every particle on the paths at every time. H2 only falls
    here (dissociation), so the range of the two ends bounds every state between them.
    The premise: at the last snapshot every particle of the domain lies inside the range
    of the reference snapshot, up to FROZEN_H2_REL.

    Raises
    ------
    GateExit
        EXIT_VOID if the n_H2 of a domain particle left the reference range.
    """
    n_first, n_last = h2_column_density(first)[mask], h2_column_density(last)[mask]
    lo_ref, hi_ref = float(n_first.min()), float(n_first.max())
    out = int(
        np.sum(
            (n_last < lo_ref * (1.0 - FROZEN_H2_REL))
            | (n_last > hi_ref * (1.0 + FROZEN_H2_REL))
        )
    )
    if out:
        raise GateExit(
            EXIT_VOID,
            f"void run, H2 changed: n_H2 of {out} probe-domain particles left the reference "
            f"range {lo_ref:.5g} to {hi_ref:.5g} cm^-3 by more than {FROZEN_H2_REL} "
            f"(last range {float(n_last.min()):.5g} to {float(n_last.max()):.5g})",
        )
    n_min, n_max = min(lo_ref, float(n_last.min())), max(hi_ref, float(n_last.max()))
    rep.info(
        "n_H2 range premise",
        f"all domain particles inside the reference range; bar range {n_min:.5g} to {n_max:.5g} cm^-3",
    )
    return n_min, n_max


def column_spread_bar(
    r: np.ndarray,
    n_h2: np.ndarray,
    mask: np.ndarray,
    temp,
    dens,
    n_range: Optional[Tuple[float, float]] = None,
) -> np.ndarray:
    """Return the bar on |ln(t / f(n_H2 r))| set by the spread of n_H2 alone.

    The carried column is a sum over the hops of the receiver's n_H2 times the path,
    and the paths telescope to r, so it lies between n_min r and n_max r (n over the
    probe domain). The held transmission is a mean of factors of such columns, and
    the factor falls with the column, so the bar is the larger of
    |ln f(n_max r) / f(n_i r)| and |ln f(n_min r) / f(n_i r)|, in the fit of the probe.
    """
    pos = mask & (n_h2 > 0.0)
    n_min, n_max = (
        n_range if n_range else (float(n_h2[pos].min()), float(n_h2[pos].max()))
    )
    base = np.log(transmission(n_h2 * r, temp, dens))
    lo = np.abs(np.log(transmission(n_max * r, temp, dens)) - base)
    hi = np.abs(np.log(transmission(n_min * r, temp, dens)) - base)
    return np.maximum(lo, hi)


def gate_starlit(args: argparse.Namespace) -> int:
    rep = Report()
    fix = load_json(args.fixture, "geometry file")
    first = load_snapshot(args.first, REQUIRED_GAS)
    last = load_snapshot(args.last, REQUIRED_GAS + REQUIRED_LW + REQUIRED_OWN)
    require_finite(
        "last snapshot",
        **{k: last[k] for k in REQUIRED_GAS + REQUIRED_LW + REQUIRED_OWN},
    )
    check_box_faces(first)
    params = read_params(args.params)
    calibrate_temperatures(first, last, params, rep)
    check_premise_parameters(params, rep, args.mode, args.startup_log)
    counters = load_json(args.counters, "counters file") if args.counters else {}

    r, h, mask = probe_geometry(last, np.array(fix["star_kpc"]))
    n_h2 = h2_column_density(last)
    column = n_h2 * r
    read = probe_premises(first, last, mask, column, rep, frozen_h2=False)
    n_range = h2_range_premise(first, last, mask, rep)
    own, meas = own_transmission(last)
    lit = read & (own > 0.0)
    rep.info(
        "probes",
        f"{int(mask.sum())} in r < L/2 - 2h, {int(read.sum())} read, {int(lit.sum())} lit",
    )
    if not lit.any():
        raise GateExit(EXIT_VOID, "void run, no lit probe")

    require_finite("measured own ratio", ratio=meas[lit])
    dens = gas_number_density(last)
    temp = last["Temperatures"]
    exp = transmission(column, temp, dens)
    slope = log_slope(np.maximum(column, 1.0), temp, dens)
    l_hop = hop_length(counters, h, rep)
    n_float = r / np.maximum(h, 1e-300) + 2.0
    float_term = n_float * (2.0 + 2.0 * slope) * FLOAT32_EPS
    bar_last = column_spread_bar(r, n_h2, mask, temp, dens) + float_term
    bar = column_spread_bar(r, n_h2, mask, temp, dens, n_range) + float_term
    widening = float(np.median(bar[lit]) / np.median(bar_last[lit]) - 1.0)
    rep.info(
        "bar widening",
        f"the bar uses the n_H2 range of the reference and the last snapshot; its median is "
        f"{100.0 * widening:.3g} per cent above the bar of the last snapshot alone",
    )
    dev = np.abs(np.log(np.maximum(meas, 1e-300) / exp))
    over = lit & (dev > bar)
    rep.info(
        "column spread",
        f"n_H2 over the probe domain from {float(n_h2[mask].min()):.4g} to "
        f"{float(n_h2[mask].max()):.4g} cm^-3; median bar {float(np.median(bar[lit])):.3g}",
    )
    rep.check(
        "own transmission against f_sh(n_H2 r)",
        not over.any(),
        f"{int(over.sum())} of {int(lit.sum())} probes outside |ln(own ratio/f)| <= the n_H2 spread bar; "
        f"worst margin {np.max((dev / bar)[lit]):.3g} of the bar; "
        f"median ln(own ratio/f) {float(np.median(np.log(np.maximum(meas[lit], 1e-300) / exp[lit]))):.4g}",
    )
    lw = last["LWSpecificEnergies"]
    mean_ratio = last["LWDissociatingSpecificEnergies"] / np.where(lw > 0, lw, 1.0)
    lit_mean = read & (lw > 0.0)
    if lit_mean.any():
        rep.reading(
            "kernel-mean ratio LWLineTransmissions",
            f"median ln(ratio/f) {float(np.median(np.log(np.maximum(mean_ratio[lit_mean], 1e-300) / exp[lit_mean]))):.4g} "
            f"over {int(lit_mean.sum())} probes (biased up near the light front)",
        )

    # Lit set: the light moves one hop per step, so the lit probes fill a ball.
    n_hops = light_hops(first, last)
    h_dom = h[mask]
    r_full = (n_hops - 1) * float(h_dom.min())
    r_reach = n_hops * KERNEL_GAMMA * float(h_dom.max())
    inside = mask & (r < r_full)
    n_full = int(inside.sum())
    n_reach = int(np.sum(mask & (r <= r_reach)))
    lit_dom = mask & (own > 0.0)
    n_lit = int(lit_dom.sum())
    rep.info(
        "light front",
        f"{n_hops} hops: every probe at r < (hops - 1) h_min = {r_full:.4g} cm ({n_full} probes) "
        f"is lit; none beyond hops gamma h_max = {r_reach:.4g} cm ({n_reach} probes)",
    )
    rep.check(
        "lit count within the light front",
        n_full <= n_lit <= n_reach,
        f"{n_lit} lit probes; bounds [{n_full}, {n_reach}]",
    )
    unlit_in = inside & ~lit_dom
    rep.check(
        "lit set complete inside the light front",
        not unlit_in.any(),
        f"{int(unlit_in.sum())} of {n_full} probes at r < {r_full:.4g} cm are unlit",
    )

    # Count gate and minimum, on the lit probes only: an unlit probe carries no transmission.
    c = int(np.sum(lit & (meas < 0.99)))
    med_t, med_n = float(np.median(temp[lit])), float(np.median(dens[lit]))
    n_h2_med = float(np.median(n_h2[lit]))
    r099 = column_for_transmission(0.99, med_t, med_n) / n_h2_med
    lh = float(np.median(l_hop[lit]))
    lo_count = int(np.sum(lit & (r >= r099 + lh)))
    hi_count = int(np.sum(lit & (r >= r099 - lh)))
    if lo_count == 0:
        raise GateExit(
            EXIT_VOID,
            "void run, no lit particle is expected to have a transmission below 0.99 (positive-control count 0)",
        )
    rep.check(
        "count of lit transmission below 0.99",
        lo_count <= c <= hi_count,
        f"{c} lit particles; shell bounds [{lo_count}, {hi_count}] around r_0.99 = {r099:.3g} cm",
    )
    deepest = np.argmax(np.where(lit, r, -1.0))
    f_deepest = float(exp[deepest])
    rep.check(
        "minimum lit transmission at or below the deepest lit probe",
        float(meas[lit].min()) <= f_deepest * math.exp(min(float(bar[deepest]), 50.0)),
        f"min {float(meas[lit].min()):.4g}, f_sh of the deepest lit probe {f_deepest:.4g} times exp(its bar)",
    )

    # N2 detection range, with the slope of the bar.
    n2 = np.abs(np.log(transmission(2.0 * column, temp, dens) / exp))
    detect = int(np.sum(lit & (n2 > bar)))
    rep.info(
        "N2 effect",
        f"|ln| of doubling n_H2 is {np.median(n2[lit]):.3g} at the median probe",
    )
    if detect == 0:
        raise GateExit(
            EXIT_VOID, "no probe can detect the doubled-column error (N2 is not a test)"
        )
    rep.info("N2 detectable probes", f"{detect} of {int(lit.sum())}")
    return rep.exit_code()


def fit_parameters(
    temp: np.ndarray, dens_cm3: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """Return the fit exponent alpha and the line width b_5 of each particle."""
    b5 = np.sqrt(2.0 * K_B * temp / (2.0 * M_H)) / 1e5
    tp = np.log10(np.clip(temp, 100.0, 8000.0))
    npr = np.log10(np.clip(dens_cm3, 1e-8, 1e7))
    alpha = (0.8711 * tp - 1.928) * np.exp(-0.2856 * npr) + (-0.9639 * tp + 3.892)
    return alpha, b5


def check_fit_uniform(snap: dict, eps: float, x_max: float, rep: Report) -> None:
    """Refuse a run whose fit (alpha, b_5) differs from cell to cell.

    First-order tolerance, needs physics review. The sign and 1/pi_l bars of the merge
    hold only where the fit is the same on every path. To first order a spread d_alpha moves ln F by at most d_alpha ln(1 + x_max / b_5)
    and a relative spread d_b of b_5 moves it by at most alpha d_b. Each must stay
    below the bar eps of the probes, or the fit is not uniform.

    Raises
    ------
    GateExit
        EXIT_VOID with 'fit not uniform'.
    """
    alpha, b5 = fit_parameters(snap["Temperatures"], gas_number_density(snap))
    d_alpha = float(alpha.max() - alpha.min())
    d_b = float((b5.max() - b5.min()) / b5.min())
    shift_a = d_alpha * math.log1p(x_max / float(b5.min()))
    shift_b = float(alpha.max()) * d_b
    rep.info(
        "fit uniformity",
        f"alpha spread {d_alpha:.3g} (ln F shift {shift_a:.3g}), "
        f"b_5 spread {d_b:.3g} (ln F shift {shift_b:.3g}), bar {eps:.3g}",
    )
    if max(shift_a, shift_b) > eps:
        raise GateExit(
            EXIT_VOID,
            "fit not uniform: the merge bars do not apply (use T = 163.4 K or a uniform gas)",
        )


def merge_rules(
    shares: List[dict], temp: float, dens: float, down: float
) -> Dict[str, float]:
    """Return the transmission after a downstream column for each merge rule.

    ``true`` keeps every share's own column and transmission, ``ii`` is the carried
    column of the code (its factor is the harmonic mean of the factors of the merged
    columns, weighted by the dissociating energy; in gas with one fit it is form (ii)),
    ``i`` the energy-weighted column and ``r3`` the transmission-weighted column. A
    share's ``t`` is its transmission at the merge; without it the share crossed gas
    of the merge fit only, t = f(N) / f(0). ``pi_l`` is the weight in the sum over the
    factor of the share with the smallest column: the bound is 1 / pi_l, which is
    1 / w_l only in gas with one fit.
    """
    w = np.array([s["w"] for s in shares], float)
    w = w / w.sum()
    n = np.array([s["N_cm2"] for s in shares], float)
    f = lambda col: float(transmission(col, temp, dens))
    t = np.array([float(s["t"]) if "t" in s else f(c) for s, c in zip(shares, n)])
    t_k = float(np.sum(w * t))
    over = w * t / np.array([f(c) for c in n])
    t_harm = t_k / float(np.sum(over))
    n_eff = column_for_transmission(t_harm, temp, dens) if t_harm < 1.0 else 0.0
    n_r3 = float(np.sum(w * t * n) / np.sum(w * t))
    lo = int(np.argmin(n))
    return {
        "true": float(sum(wi * ti * f(c + down) / f(c) for wi, ti, c in zip(w, t, n))),
        "ii": t_k * f(n_eff + down) / f(n_eff),
        "i": f(float(np.sum(w * n)) + down),
        "r3": t_k * f(n_r3 + down) / f(n_r3),
        "n_eff": n_eff,
        "w_l": float(w[lo]),
        "pi_l": float(over[lo] / np.sum(over)),
    }


TRUE_TRANSMISSION_READ_MIN = 1e-2


def gate_twosource(args: argparse.Namespace) -> int:
    rep = Report()
    fix = load_json(args.fixture, "geometry file")
    first = load_snapshot(args.first, REQUIRED_GAS)
    last = load_snapshot(args.last, REQUIRED_GAS + REQUIRED_LW + REQUIRED_OWN)
    require_finite(
        "last snapshot",
        **{k: last[k] for k in REQUIRED_GAS + REQUIRED_LW + REQUIRED_OWN},
    )
    check_box_faces(first)
    params = read_params(args.params)
    calibrate_temperatures(first, last, params, rep)
    check_premise_parameters(params, rep, args.mode, args.startup_log)
    counters = load_json(args.counters, "counters file") if args.counters else {}
    rep.info(
        "held light at the merge",
        "0 in every run (the end of the force loop zeroes the held energy), "
        "so the held term of the merge is not exercised here",
    )
    ids = {int(i): k for k, i in enumerate(last["ParticleIDs"])}
    probes = [p for p in fix["probes"] if int(p["id"]) in ids]
    idx = np.array([ids[int(p["id"])] for p in probes], dtype=int)
    mask = np.zeros(len(last["ParticleIDs"]), bool)
    mask[idx] = True
    read = probe_premises(first, last, mask, np.full(len(mask), np.inf), rep)
    dens, temp = gas_number_density(last), last["Temperatures"]
    lw, meas = own_transmission(last)
    require_finite("measured own ratio", ratio=meas[idx])
    cm = last["unit_length"]
    l_hop_all = hop_length(counters, last["SmoothingLengths"] * cm, rep)

    rows = []
    for p, k in zip(probes, idx):
        if not read[k] or lw[k] <= 0.0:
            continue
        down = float(p.get("downstream_N_cm2", 0.0))
        rules = merge_rules(p["shares"], temp[k], dens[k], down)
        n_probe = rules["n_eff"] + down
        s = float(log_slope(max(n_probe, 1.0), temp[k], dens[k]))
        eps_hop = s * l_hop_all[k] / max(float(p["r_b_kpc"]) * cm, l_hop_all[k])
        eps_float = float(p.get("n_merges", 1)) * (2.0 + 2.0 * s) * FLOAT32_EPS
        eps = eps_hop + eps_float
        r_meas = meas[k] / rules["true"]
        r_pred = rules["ii"] / rules["true"]
        r3_fails = rules["r3"] / rules["true"] > (1.0 / rules["pi_l"]) * (1.0 + eps)
        i_fails = rules["i"] / rules["true"] < 1.0 - eps
        rows.append(
            (
                p["id"],
                r_meas,
                r_pred,
                rules["pi_l"],
                eps,
                r3_fails,
                i_fails,
                rules["true"],
            )
        )
    if not rows:
        raise GateExit(EXIT_VOID, "void run, no lit merge probe")
    a = np.array([[r[1], r[2], r[3], r[4]] for r in rows])
    r_meas, r_pred, pi_l, eps = a.T
    t_true = np.array([r[7] for r in rows])
    check_fit_uniform(first, float(np.min(eps)), 1e24 / 5e14, rep)
    main_bad = np.abs(r_meas / r_pred - 1.0) > eps
    sign_bad = r_meas < 1.0 - eps
    bound_bad = r_meas > (1.0 / pi_l) * (1.0 + eps)
    rep.check(
        "ratio equals the form (ii) prediction",
        not main_bad.any(),
        f"{int(main_bad.sum())} of {len(rows)} probes outside |R/R_pred - 1| <= eps_hop + n_m eps_m; "
        f"worst margin {np.max(np.abs(r_meas / r_pred - 1.0) / eps):.3g} of the bar",
    )
    # Readings only: the sign and the bound hold where the exact transmission is
    # above TRUE_TRANSMISSION_READ_MIN; they never change the exit code.
    seen = t_true > TRUE_TRANSMISSION_READ_MIN
    n_seen = int(seen.sum())
    if n_seen:
        rep.info(
            "sign reading: R >= 1 - eps",
            f"{int((sign_bad & seen).sum())} of {n_seen} probes with exact transmission "
            f"above {TRUE_TRANSMISSION_READ_MIN:g} below; smallest R {r_meas[seen].min():.6g}",
        )
        rep.info(
            "bound reading: R <= (1 / pi_l)(1 + eps)",
            f"{int((bound_bad & seen).sum())} of {n_seen} probes above; "
            f"largest R {r_meas[seen].max():.6g}",
        )
    else:
        rep.info(
            "sign and bound readings",
            f"no probe has an exact transmission above {TRUE_TRANSMISSION_READ_MIN:g}",
        )
    n_bound = sum(1 for r in rows if r[5])
    n_sign = sum(1 for r in rows if r[6])
    rep.info(
        "discrimination",
        f"{n_bound} of {len(rows)} probes would fail the bound for the r3 rule, "
        f"{n_sign} would fail the sign bar for form (i); the others are too shallow to discriminate",
    )
    if n_bound == 0:
        raise GateExit(
            EXIT_VOID,
            "no probe is deep enough to discriminate the bound (downstream column too small)",
        )
    return rep.exit_code()


def velocity_bar(
    r: np.ndarray,
    h: np.ndarray,
    depth: np.ndarray,
    n_range: Tuple[float, float],
    reach: float,
    temp: np.ndarray,
    dens: np.ndarray,
    reset: bool,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return the lower and upper bound of ln(own ratio) of the cloud probes.

    With a reset the carried column is n times the path from the sender of the resetting
    hop, which lies up to one reach in front of the face: n_min d <= N <= n_max (d + R),
    with d the path inside the cloud and R = gamma h_max the largest hop. Without a reset
    the column telescopes to n r over all the gas: n_min r <= N <= n_max r. The factor
    falls with the column, so the bound on ln t is that of the opposite column. The float
    term is that of the star-lit gate for a slope of at most 1.
    """
    n_min, n_max = n_range
    eps = (r / np.maximum(h, 1e-300) + 2.0) * 4.0 * FLOAT32_EPS
    col_hi = n_max * (depth + reach) if reset else n_max * r
    col_lo = n_min * depth if reset else n_min * r
    lo = np.log(transmission(col_hi, temp, dens)) - eps
    hi = np.log(transmission(col_lo, temp, dens)) + eps
    return lo, hi


def gate_velocity(args: argparse.Namespace) -> int:
    rep = Report()
    fix = load_json(args.fixture, "geometry file")
    first = load_snapshot(args.first, REQUIRED_GAS)
    last = load_snapshot(args.last, REQUIRED_GAS + REQUIRED_LW + REQUIRED_OWN)
    require_finite(
        "last snapshot",
        **{k: last[k] for k in REQUIRED_GAS + REQUIRED_LW + REQUIRED_OWN},
    )
    check_box_faces(first)
    params = read_params(args.params)
    calibrate_temperatures(first, last, params, rep)
    check_premise_parameters(params, rep, args.mode, args.startup_log)
    star = np.array(fix["star_kpc"])
    r, h, mask = probe_geometry(last, star)
    cm = last["unit_length"]
    face, outer = (v * cm for v in fix["cloud_radius_kpc"])
    in_cloud = mask & (r > face) & (r < outer)
    read = probe_premises(first, last, in_cloud, np.full(len(mask), np.inf), rep)
    lw, meas = own_transmission(last)
    lit = read & (lw > 0.0)
    if not lit.any():
        raise GateExit(EXIT_VOID, "void run, no lit probe inside the cloud")
    require_finite("measured own ratio", ratio=meas[lit])
    dens, temp = gas_number_density(last), last["Temperatures"]
    n_h2 = h2_column_density(last)
    require_finite("n_H2", n_h2=n_h2[mask])
    dv = float(fix["cloud_speed_widths"])
    reset = dv > WIDTH_RESET if args.case == "auto" else args.case == "reset"
    gas = in_cloud if reset else mask
    n_range = (float(n_h2[gas].min()), float(n_h2[gas].max()))
    near = mask & (r > face - 2.0 * KERNEL_GAMMA * float(h[mask].max()))
    reach = KERNEL_GAMMA * float(h[near].max())
    depth = np.maximum(r - face, 0.0)
    lo, hi = velocity_bar(r, h, depth, n_range, reach, temp, dens, reset)
    ln_t = np.log(np.maximum(meas, 1e-300))
    require_finite("velocity bar", lo=lo[lit], hi=hi[lit], ln_ratio=ln_t[lit])
    dark = lit & (ln_t < lo)
    bright = lit & (ln_t > hi)
    shallow = lit & (depth < h)
    rep.info(
        "velocity bar",
        f"n_H2 {n_range[0]:.5g} to {n_range[1]:.5g} cm^-3, reach gamma h_max = {reach:.4g} cm, "
        f"{int(lit.sum())} lit probes, {int(shallow.sum())} less than one h deep",
    )
    rep.check(
        f"cloud probes follow the {'own' if reset else 'full'} column",
        not (dark | bright).any(),
        f"{int(dark.sum())} too dark and {int(bright.sum())} too bright of {int(lit.sum())} probes "
        f"({int((dark | bright)[shallow].sum())} of {int(shallow.sum())} less than one h deep); "
        f"cloud speed {dv:g} widths, expected {'reset' if reset else 'no reset'}",
    )
    return rep.exit_code()


def raw_gas_datasets(path: str) -> Dict[str, np.ndarray]:
    """Return every PartType0 dataset of a snapshot as stored, sorted by particle ID."""
    try:
        f = h5py.File(path, "r")
    except OSError as err:
        raise GateExit(EXIT_MISSING, f"snapshot {path} cannot be opened: {err}")
    with f:
        if "PartType0" not in f or "ParticleIDs" not in f["PartType0"]:
            raise GateExit(EXIT_MISSING, f"{path} has no PartType0 particle IDs")
        g = f["PartType0"]
        order = np.argsort(g["ParticleIDs"][:])
        return {k: g[k][:][order] for k in g}


def differing_rows(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Return a mask of the particles whose stored bytes differ."""
    bx = np.ascontiguousarray(x).view(np.uint8).reshape(len(x), -1)
    by = np.ascontiguousarray(y).view(np.uint8).reshape(len(y), -1)
    return np.any(bx != by, axis=1)


def kernel_neighbours(snap: dict) -> List[List[int]]:
    """Return, for each particle, the particles within gamma h of it (periodic)."""
    from scipy.spatial import cKDTree

    box = snap["box"]
    pos = np.mod(snap["Coordinates"], box)
    tree = cKDTree(pos, boxsize=box)
    return tree.query_ball_point(pos, KERNEL_GAMMA * snap["SmoothingLengths"])


def gate_parity(args: argparse.Namespace) -> int:
    """Compare two runs particle by particle: two rank counts, or a restart.

    restart: the same binary and rank count, resumed from a restart file. Every
    LW and packet field must be equal bit for bit.

    ranks: the reordered float sums of the transport give a round-off gap. A
    particle's field is a sum of positive shares, each a weight times the
    sender's energy. A sum of k float32 terms errs by at most k 2^-24 of its
    value, and a smoothing length that differs by delta_h changes a weight of
    power p by p delta_h of its value away from the edge of the kernel. Over n
    steps the light crosses n such sums, so a particle that held light at the
    previous snapshot has a gap of at most n (k 2^-24 + p delta_h) of its own
    value. A particle first lit after the previous snapshot is at the light
    front: its light came through weights near the kernel edge, where the
    relative change of a weight is not bounded, so it is held to the same bar
    times the largest value in its kernel, which is weaker.
    """
    rep = Report()
    gated = ("LWSpecificEnergies", "LWDissociatingSpecificEnergies", PARITY_COLUMN)
    geom = ("ParticleIDs", "Coordinates", "SmoothingLengths")
    a = load_snapshot(args.a, geom + gated)
    b = load_snapshot(args.b, geom + gated)
    counters = load_json(args.counters, "counters file")
    need = (
        ["probes_on_other_rank", "foreign_bookings", "steps"]
        if args.mode == "ranks"
        else ["held_lw_slots_at_restart", "refits_after_restart"]
    )
    for key in need:
        if not float(counters.get(key, 0)) > 0:
            raise GateExit(
                EXIT_VOID,
                f"void run, {key} is {counters.get(key, 'absent')} (must be above 0)",
            )
    if not a["time"] == b["time"]:
        raise GateExit(
            EXIT_VOID,
            f"void run, the snapshots are at times {a['time']} and {b['time']}",
        )
    oa, ob = np.argsort(a["ParticleIDs"]), np.argsort(b["ParticleIDs"])
    if not np.array_equal(a["ParticleIDs"][oa], b["ParticleIDs"][ob]):
        raise GateExit(
            EXIT_VOID, "void run, the two snapshots hold different particles"
        )
    for key in gated:
        require_finite(key, a=a[key], b=b[key])
    rep.info("premises", ", ".join(f"{k} {counters[k]}" for k in need))

    if args.mode == "restart":
        if args.stop is None:
            raise GateExit(EXIT_MISSING, "restart mode needs --stop")
        stop = load_snapshot(args.stop, ("ParticleIDs", "LWSpecificEnergies"))
        os_ = np.argsort(stop["ParticleIDs"])
        if not np.array_equal(stop["ParticleIDs"][os_], a["ParticleIDs"][oa]):
            raise GateExit(
                EXIT_VOID, "void run, the stop snapshot holds other particles"
            )
        dark = stop["LWSpecificEnergies"][os_] == 0.0
        new_lit = int(np.sum(dark & (a["LWSpecificEnergies"][oa] > 0.0)))
        if new_lit == 0:
            raise GateExit(
                EXIT_VOID, "void run, no particle dark at the stop is lit after it"
            )
        rep.info("lit after the restart", f"{new_lit} particles dark at the stop")
        ra, rb = raw_gas_datasets(args.a), raw_gas_datasets(args.b)
        for key in sorted(set(ra) | set(rb)):
            if key not in ra or key not in rb:
                raise GateExit(EXIT_MISSING, f"{key} is in one snapshot only")
            neq = differing_rows(ra[key], rb[key])
            detail = f"{int(neq.sum())} of {len(neq)} particles differ"
            if key in gated or key.startswith(("LW", "Packet")):
                rep.check(f"{key} bitwise", not neq.any(), detail)
            else:
                rep.reading(f"{key} bitwise", detail)
        return rep.exit_code()

    steps = float(counters["steps"])
    ngb = kernel_neighbours(a)
    k_terms = max(len(n) for n in ngb)
    params = read_params(args.params) if args.params else {}
    p_power = float(params.get("GEARFeedback:ISRF_packet_cone_power", 4.0)) + 4.0
    ha, hb = a["SmoothingLengths"][oa], b["SmoothingLengths"][ob]
    delta_h = float(np.max(np.abs(ha - hb) / np.maximum(ha, hb)))
    bar = steps * (k_terms * FLOAT32_EPS + p_power * delta_h)
    rep.info(
        "round-off bar",
        f"n (k 2^-24 + p delta_h) = {steps:g} ({k_terms} x 2^-24 + {p_power:g} x "
        f"{delta_h:.3g}) = {bar:.3g} of the largest value in the kernel",
    )
    if args.previous is None:
        raise GateExit(EXIT_MISSING, "ranks mode needs --previous")
    prev = load_snapshot(args.previous, ("ParticleIDs", "LWSpecificEnergies"))
    op = np.argsort(prev["ParticleIDs"])
    if not np.array_equal(prev["ParticleIDs"][op], a["ParticleIDs"][oa]):
        raise GateExit(
            EXIT_VOID, "void run, the previous snapshot holds other particles"
        )
    held = np.empty(len(oa), bool)
    held[oa] = prev["LWSpecificEnergies"][op] > 0.0
    front = ~held & (a["LWSpecificEnergies"] > 0.0)
    if not held.any():
        raise GateExit(EXIT_VOID, "void run, no particle held light at --previous")
    rep.info(
        "particles",
        f"{int(held.sum())} held light at the previous snapshot, {int(front.sum())} "
        "first lit after it (light front)",
    )
    for key in gated:
        x = a[key]
        y = np.empty_like(x)
        y[oa] = b[key][ob]
        scale = np.maximum(np.abs(x), np.abs(y))
        gap = np.abs(x - y)
        own = np.where(scale > 0, gap / np.where(scale > 0, scale, 1), 0.0)
        over_held = held & (gap > bar * scale)
        rep.check(
            f"{key} round-off, lit before",
            not over_held.any(),
            f"{int(over_held.sum())} of {int(held.sum())} particles over the bar, "
            f"largest gap over the own value {own[held].max():.3g}, bar {bar:.3g}",
        )
        local = np.array([scale[n].max() for n in ngb])
        rest = ~held
        over_rest = rest & (gap > bar * local)
        worst = np.where(local > 0, gap / np.where(local > 0, local, 1), 0.0)
        rep.check(
            f"{key} round-off, light front and dark",
            not over_rest.any(),
            f"{int(over_rest.sum())} of {int(rest.sum())} particles over the bar times "
            f"the kernel maximum, largest {worst[rest].max():.3g}; largest gap over the "
            f"own value {own[rest].max():.3g} (read only)",
        )
    return rep.exit_code()


def gate_refusal(args: argparse.Namespace) -> int:
    rep = Report()
    try:
        text = open(args.log).read()
    except OSError as err:
        raise GateExit(EXIT_MISSING, f"log {args.log} is missing: {err}")
    pattern = (
        r"^\[\d+\.\d+\] \S+:\S+\(\):\d+: GEARFeedback:ISRF_H2_line_absorption: 1 is not supported"
        if args.expect == "shielded_sum_refused"
        else args.pattern
    )
    found = re.search(pattern, text, re.M) is not None
    rep.check(
        "run stopped at start-up",
        args.exit_status != 0,
        f"exit status {args.exit_status}",
    )
    rep.check("refusal text in the log", found, "anchored on SWIFT's error() format")
    return rep.exit_code()


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp: argparse.ArgumentParser) -> None:
        sp.add_argument("--first", required=True, help="First snapshot")
        sp.add_argument("--last", required=True, help="Last snapshot")
        sp.add_argument("--fixture", default="fixture.json")
        sp.add_argument("--params", default="used_parameters.yml")
        sp.add_argument(
            "--mode",
            choices=["a", "c"],
            default="a",
            help="a: key in used_parameters.yml; c: start-up log",
        )
        sp.add_argument("--startup-log", default=None)
        sp.add_argument("--counters", default=None)

    for name, fn in (
        ("starlit", gate_starlit),
        ("twosource", gate_twosource),
        ("velocity", gate_velocity),
    ):
        sp = sub.add_parser(name)
        common(sp)
        if name == "velocity":
            sp.add_argument(
                "--case", choices=["auto", "reset", "noreset"], default="auto"
            )
        sp.set_defaults(fn=fn)
    sp = sub.add_parser("parity")
    sp.add_argument("--a", required=True)
    sp.add_argument("--b", required=True)
    sp.add_argument("--mode", choices=["ranks", "restart"], required=True)
    sp.add_argument("--counters", required=True)
    sp.add_argument(
        "--stop", default=None, help="restart: snapshot of run a at the stop time"
    )
    sp.add_argument(
        "--params", default=None, help="ranks: used_parameters.yml, for the cone power"
    )
    sp.add_argument(
        "--previous",
        default=None,
        help="ranks: snapshot of run a before the compared one, to find the light front",
    )
    sp.set_defaults(fn=gate_parity)
    sp = sub.add_parser("refusal")
    sp.add_argument("--log", required=True)
    sp.add_argument("--exit-status", type=int, required=True)
    sp.add_argument("--expect", default="shielded_sum_refused")
    sp.add_argument("--pattern", default="")
    sp.set_defaults(fn=gate_refusal)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    try:
        code = args.fn(args)
    except GateExit as stop:
        print(f"FAIL  {stop.message}")
        code = stop.code
    print(f"exit class {code}")
    return code


if __name__ == "__main__":
    sys.exit(main())
