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

REQUIRED_LW = (
    "LWSpecificEnergies",
    "LWDissociatingSpecificEnergies",
    "LWLineTransmissions",
)
REQUIRED_GAS = (
    "ParticleIDs",
    "Coordinates",
    "SmoothingLengths",
    "Densities",
    "Temperatures",
    "H2I",
)


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


def load_json(path: str, what: str) -> dict:
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
        out["unit_length"] = float(
            np.asarray(f["Units"].attrs["Unit length in cgs (U_L)"]).ravel()[0]
        )
        out["unit_mass"] = float(
            np.asarray(f["Units"].attrs["Unit mass in cgs (U_M)"]).ravel()[0]
        )
    return out


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
    first: dict, last: dict, mask: np.ndarray, column: np.ndarray, rep: Report
) -> np.ndarray:
    """Check frozen H2, the LW source and the temperature; return the read mask.

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
        for k in ("H2I", "Temperatures"):
            tmp = np.empty_like(first[k])
            tmp[order1] = first[k][order0]
            first[k] = tmp
    dx = np.abs(last["H2I"] - first["H2I"]) / np.maximum(first["H2I"], 1e-300)
    bad_h2 = int(np.sum(mask & (dx > FROZEN_H2_REL)))
    if bad_h2:
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


def gate_starlit(args: argparse.Namespace) -> int:
    rep = Report()
    fix = load_json(args.fixture, "geometry file")
    first = load_snapshot(args.first, REQUIRED_GAS)
    last = load_snapshot(args.last, REQUIRED_GAS + REQUIRED_LW)
    require_finite("last snapshot", **{k: last[k] for k in REQUIRED_GAS + REQUIRED_LW})
    params = read_params(args.params)
    check_premise_parameters(params, rep, args.mode, args.startup_log)
    counters = load_json(args.counters, "counters file") if args.counters else {}

    r, h, mask = probe_geometry(last, np.array(fix["star_kpc"]))
    n_h2 = h2_column_density(last)
    column = n_h2 * r
    read = probe_premises(first, last, mask, column, rep)
    lw = last["LWSpecificEnergies"]
    lit = read & (lw > 0.0)
    rep.info(
        "probes",
        f"{int(mask.sum())} in r < L/2 - 2h, {int(read.sum())} read, {int(lit.sum())} lit",
    )
    if not lit.any():
        raise GateExit(EXIT_VOID, "void run, no lit probe")

    meas = last["LWDissociatingSpecificEnergies"] / np.where(lw > 0, lw, 1.0)
    require_finite("measured ratio", ratio=meas[lit])
    dens = gas_number_density(last)
    temp = last["Temperatures"]
    exp = transmission(column, temp, dens)
    slope = log_slope(np.maximum(column, 1.0), temp, dens)
    l_hop = hop_length(counters, h, rep)
    bar = slope * l_hop / r
    dev = np.abs(np.log(np.maximum(meas, 1e-300) / exp))
    over = lit & (dev > bar)
    rep.check(
        "transmission against f_sh(n_H2 r)",
        not over.any(),
        f"{int(over.sum())} of {int(lit.sum())} probes outside |ln(ratio/f)| <= slope l_hop / r; "
        f"worst margin {np.max((dev / bar)[lit]):.3g} of the bar",
    )

    # Count gate.
    t_line = last["LWLineTransmissions"]
    in_dom = mask
    c = int(np.sum(in_dom & (t_line < 0.99)))
    med_t, med_n = float(np.median(temp[mask])), float(np.median(dens[mask]))
    n_h2_med = float(np.median(n_h2[mask]))
    r099 = column_for_transmission(0.99, med_t, med_n) / n_h2_med
    lh = float(np.median(l_hop[mask]))
    lo_count = int(np.sum(in_dom & (r >= r099 + lh)))
    hi_count = int(np.sum(in_dom & (r >= r099 - lh)))
    if c == 0:
        raise GateExit(
            EXIT_VOID,
            "void run, no particle has a transmission below 0.99 (positive-control count 0)",
        )
    rep.check(
        "count of transmission below 0.99",
        lo_count <= c <= hi_count,
        f"{c} particles; shell bounds [{lo_count}, {hi_count}] around r_0.99 = {r099:.3g} cm",
    )
    deepest = np.argmax(np.where(mask, r, -1.0))
    f_deepest = float(exp[deepest])
    rep.check(
        "minimum transmission at or below the deepest probe",
        float(t_line[in_dom].min())
        <= f_deepest * math.exp(min(float(bar[deepest]), 50.0)),
        f"min {float(t_line[in_dom].min()):.4g}, f_sh of the deepest probe {f_deepest:.4g} times exp(its hop bar)",
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

    The sign and 1/w0 bars of the merge hold only where the fit is the same on every
    path. To first order a spread d_alpha moves ln F by at most d_alpha ln(1 + x_max / b_5)
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

    ``true`` keeps every share's own column, ``ii`` stores the exact mean transmission,
    ``i`` the energy-weighted column and ``r3`` the transmission-weighted column.
    """
    w = np.array([s["w"] for s in shares], float)
    w = w / w.sum()
    n = np.array([s["N_cm2"] for s in shares], float)
    f = lambda col: float(transmission(col, temp, dens))
    t = np.array([f(c) for c in n])
    t_k = float(np.sum(w * t))
    n_eff = column_for_transmission(t_k, temp, dens) if t_k < 1.0 else 0.0
    n_r3 = float(np.sum(w * t * n) / np.sum(w * t))
    return {
        "true": float(sum(wi * f(c + down) for wi, c in zip(w, n))),
        "ii": f(n_eff + down),
        "i": f(float(np.sum(w * n)) + down),
        "r3": t_k * f(n_r3 + down) / f(n_r3),
        "n_eff": n_eff,
        "w_l": float(w[np.argmin(n)]),
    }


def gate_twosource(args: argparse.Namespace) -> int:
    rep = Report()
    fix = load_json(args.fixture, "geometry file")
    first = load_snapshot(args.first, REQUIRED_GAS)
    last = load_snapshot(args.last, REQUIRED_GAS + REQUIRED_LW)
    require_finite("last snapshot", **{k: last[k] for k in REQUIRED_GAS + REQUIRED_LW})
    params = read_params(args.params)
    check_premise_parameters(params, rep, args.mode, args.startup_log)
    counters = load_json(args.counters, "counters file")
    if float(counters.get("held_lw_slots_at_merge", 0)) <= 0:
        raise GateExit(
            EXIT_VOID,
            "void run, no slot holds LW energy at the merge (held_lw_slots_at_merge is 0)",
        )
    ids = {int(i): k for k, i in enumerate(last["ParticleIDs"])}
    probes = [p for p in fix["probes"] if int(p["id"]) in ids]
    idx = np.array([ids[int(p["id"])] for p in probes], dtype=int)
    mask = np.zeros(len(last["ParticleIDs"]), bool)
    mask[idx] = True
    read = probe_premises(first, last, mask, np.full(len(mask), np.inf), rep)
    dens, temp = gas_number_density(last), last["Temperatures"]
    lw = last["LWSpecificEnergies"]
    meas = last["LWDissociatingSpecificEnergies"] / np.where(lw > 0, lw, 1.0)
    require_finite("measured ratio", ratio=meas[idx])
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
        r3_fails = rules["r3"] / rules["true"] > (1.0 / rules["w_l"]) * (1.0 + eps)
        i_fails = rules["i"] / rules["true"] < 1.0 - eps
        rows.append((p["id"], r_meas, r_pred, rules["w_l"], eps, r3_fails, i_fails))
    if not rows:
        raise GateExit(EXIT_VOID, "void run, no lit merge probe")
    a = np.array([[r[1], r[2], r[3], r[4]] for r in rows])
    r_meas, r_pred, w_l, eps = a.T
    check_fit_uniform(first, float(np.min(eps)), 1e24 / 5e14, rep)
    main_bad = np.abs(r_meas / r_pred - 1.0) > eps
    sign_bad = r_meas < 1.0 - eps
    bound_bad = r_meas > (1.0 / w_l) * (1.0 + eps)
    rep.check(
        "ratio equals the form (ii) prediction",
        not main_bad.any(),
        f"{int(main_bad.sum())} of {len(rows)} probes outside |R/R_pred - 1| <= eps_hop + n_m eps_m; "
        f"worst margin {np.max(np.abs(r_meas / r_pred - 1.0) / eps):.3g} of the bar",
    )
    rep.check(
        "sign: R >= 1 - eps",
        not sign_bad.any(),
        f"{int(sign_bad.sum())} of {len(rows)} probes below; smallest R {r_meas.min():.6g}",
    )
    rep.check(
        "bound: R <= (1 / w_l)(1 + eps)",
        not bound_bad.any(),
        f"{int(bound_bad.sum())} of {len(rows)} probes above; largest R {r_meas.max():.6g}",
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


def path_through(a: float, b: float, x_star: float, x_p: float, r: float) -> float:
    """Return the length of the ray star to probe inside the slab a < x < b (along x)."""
    dx = x_p - x_star
    if dx <= 0:
        return 0.0
    lo, hi = max(a, x_star), min(b, x_p)
    return r * max(0.0, hi - lo) / dx


def gate_velocity(args: argparse.Namespace) -> int:
    rep = Report()
    fix = load_json(args.fixture, "geometry file")
    first = load_snapshot(args.first, REQUIRED_GAS)
    last = load_snapshot(args.last, REQUIRED_GAS + REQUIRED_LW)
    require_finite("last snapshot", **{k: last[k] for k in REQUIRED_GAS + REQUIRED_LW})
    params = read_params(args.params)
    check_premise_parameters(params, rep, args.mode, args.startup_log)
    counters = load_json(args.counters, "counters file") if args.counters else {}
    star = np.array(fix["star_kpc"])
    r, h, mask = probe_geometry(last, star)
    cm = last["unit_length"]
    d_min = last["Coordinates"] - star[None, :]
    d_min -= last["box"] * np.round(d_min / last["box"])
    x_rel = d_min[:, 0] * cm
    slab = [(v - star[0]) * cm for v in fix["slab_x_kpc"]]
    cloud = [(v - star[0]) * cm for v in fix["cloud_x_kpc"]]
    in_cloud = mask & (x_rel > cloud[0]) & (x_rel < cloud[1])
    read = probe_premises(first, last, in_cloud, np.full(len(mask), np.inf), rep)
    lw = last["LWSpecificEnergies"]
    lit = read & (lw > 0.0)
    if not lit.any():
        raise GateExit(EXIT_VOID, "void run, no lit probe inside the cloud")
    meas = last["LWDissociatingSpecificEnergies"] / np.where(lw > 0, lw, 1.0)
    require_finite("measured ratio", ratio=meas[lit])
    dens, temp = gas_number_density(last), last["Temperatures"]
    n_h2 = float(fix["n_h2_cm3"])
    l_hop = hop_length(counters, h, rep)
    dv = float(fix["cloud_speed_widths"])
    reset = dv > WIDTH_RESET if args.case == "auto" else args.case == "reset"
    col = np.zeros(len(r))
    bar = np.zeros(len(r))
    for k in np.where(lit)[0]:
        in_cloud_path = path_through(cloud[0], cloud[1], 0.0, x_rel[k], r[k])
        in_slab_path = path_through(slab[0], slab[1], 0.0, x_rel[k], r[k])
        col[k] = n_h2 * (in_cloud_path if reset else in_cloud_path + in_slab_path)
        bar[k] = (
            float(log_slope(max(col[k], 1.0), temp[k], dens[k]))
            * l_hop[k]
            / max(in_cloud_path, l_hop[k])
        )
    exp = transmission(col, temp, dens)
    dev = np.abs(np.log(np.maximum(meas, 1e-300) / exp))
    over = lit & (dev > bar)
    rep.check(
        f"cloud probes follow the {'own' if reset else 'full'} column",
        not over.any(),
        f"{int(over.sum())} of {int(lit.sum())} probes outside the bar; "
        f"cloud speed {dv:g} widths, expected {'reset' if reset else 'no reset'}",
    )
    return rep.exit_code()


def gate_parity(args: argparse.Namespace) -> int:
    rep = Report()
    a = load_snapshot(
        args.a, ("ParticleIDs", "LWSpecificEnergies", "LWDissociatingSpecificEnergies")
    )
    b = load_snapshot(
        args.b, ("ParticleIDs", "LWSpecificEnergies", "LWDissociatingSpecificEnergies")
    )
    counters = load_json(args.counters, "counters file")
    need = (
        ["probes_on_other_rank", "foreign_bookings"]
        if args.mode == "ranks"
        else ["held_lw_slots_at_restart"]
    )
    for key in need:
        if float(counters.get(key, 0)) <= 0:
            raise GateExit(
                EXIT_VOID,
                f"void run, {key} is {counters.get(key, 'absent')} (must be above 0)",
            )
    oa, ob = np.argsort(a["ParticleIDs"]), np.argsort(b["ParticleIDs"])
    if not np.array_equal(a["ParticleIDs"][oa], b["ParticleIDs"][ob]):
        raise GateExit(
            EXIT_VOID, "void run, the two snapshots hold different particles"
        )
    for key in ("LWSpecificEnergies", "LWDissociatingSpecificEnergies"):
        require_finite(key, a=a[key], b=b[key])
        x, y = a[key][oa], b[key][ob]
        if args.mode == "restart":
            neq = x.view(np.uint64) != y.view(np.uint64)
            rep.check(
                f"{key} bitwise",
                not neq.any(),
                f"{int(neq.sum())} of {len(x)} particles differ",
            )
        else:
            k_terms = float(counters.get("max_terms_per_particle", 0))
            if k_terms <= 0:
                raise GateExit(
                    EXIT_VOID, "void run, max_terms_per_particle is absent or 0"
                )
            scale = np.maximum(np.abs(x), np.abs(y))
            rel = np.abs(x - y) / np.where(scale > 0, scale, 1.0)
            bar = k_terms * FLOAT32_EPS
            rep.info(
                f"{key} bitwise",
                f"{int(np.sum(x.view(np.uint64) != y.view(np.uint64)))} of {len(x)} particles differ",
            )
            rep.check(
                f"{key} round-off",
                bool(np.all(rel <= bar)),
                f"worst relative difference {rel.max():.3g}, bar {bar:.3g} = {k_terms:g} terms x 2^-24",
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
