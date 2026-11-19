#!/usr/bin/env python3
"""Check the H2 line absorption of the ISRF transport on the three H2 fixtures.

Subcommands
-----------
starlit     Transmission at depth against the shielding function (star-lit box).
twosource   Two stars whose lights never share a slot: the two-star run is the sum of the
            one-star runs, and each star's column follows the straight ray.
velocity    Reset of the carried column by a velocity shift.
parity      Per-particle comparison of two runs (two ranks, or restart).
stale_key   The removed key GeajayCooling:h2_transported_lw_rate_shielding is ignored.

Exit codes
----------
0 all checks passed
1 a bar failed
2 a required file, field or key is missing
3 a value is not finite
4 the run is void: a premise or a positive control is not met
5 parity only: the comparison passed without a counters file, so the premises of the
  run (foreign bookings, refits) were not checked

A missing field, a void run and a non-finite value never pass and never crash.

The shielding function is the fit of the design: x = N_H2 / 5e14 cm^-2,
b_5 = sqrt(2 k T / m_H2) / (1 km/s), f = 0.965 / (1 + x / b_5)^alpha +
0.035 exp(-8.5e-4 sqrt(1 + x)) / sqrt(1 + x), with alpha from the temperature
and density fit of Wolcott-Green and Haiman (2019). The density is floored
at 1e-8 cm^-3 as in the module. The expected transmission is f(N) / f(0).
"""

import argparse
import itertools
import json
import math
import os
import re
import sys
from typing import Dict, Iterator, List, Optional, Tuple

import h5py
import numpy as np

EXIT_OK, EXIT_BAR, EXIT_MISSING, EXIT_NONFINITE, EXIT_VOID, EXIT_REDUCED = (
    0,
    1,
    2,
    3,
    4,
    5,
)

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
REQUIRED_COLUMN = (PARITY_COLUMN,)
REQUIRED_KERNEL = ("Masses", "Velocities")
FLOAT32_TINY = float(np.finfo(np.float32).tiny)
KERNEL_WEIGHT_AT_ZERO = 21.0 / (2.0 * math.pi * KERNEL_GAMMA**3)
# Piecewise polynomial of src/kernel_hydro.h for Wendland C2 in 3D: W = Horner(c, t) with t = u - origin
# on the interval index = int(u * 4 / gamma_code); the last interval ends at the support gamma_code.
KERNEL_CODE_GAMMA = 1.9364917278289795
KERNEL_POLY_ORIGIN = np.array(
    [
        0.0,
        0.9682458639144897,
        1.4523687362670898,
        1.9364917278289795,
        1.9364917278289795,
    ]
)
KERNEL_POLY_COEFFS = np.array(
    [
        [
            0.0676042888394359,
            -0.4909317978874823,
            1.2675804873830907,
            -1.2273295640873905,
            0.0,
            0.4602486125460297,
        ],
        [
            0.0676042888394359,
            -0.16364393262916077,
            0.0,
            0.3068323910218476,
            -0.2970891935218974,
            0.08629661485238058,
        ],
        [
            0.0676042888394359,
            -2.014764810783741e-08,
            -0.15844756092288392,
            0.15341622384355558,
            -0.055704242073993125,
            0.007191387891262731,
        ],
        [0.0676042888394359, 0.16364393262916077, 0.0, 0.0, 0.0, 0.0],
        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
    ]
)
# Relative error of u = r / h in float, in units of eps32, after the float positions dx: the
# subtraction of the two float positions (1 eps in r), r^2 as three products and two sums (3 eps in
# r^2, so 1.5 in r), sqrtf (1 eps), the product with the reciprocal of h under -freciprocal-math
# (2 eps): 1 + 1.5 + 1 + 2, rounded up.
ARGUMENT_ROUNDING_EPS = 6.0
# Absolute error of u from the float positions, in units of eps32 C / h. The force loops cast
# (x - loc) to float for each particle, with loc the corner of the cell of pj (shifted for a periodic
# pair), so the cast error of a coordinate is eps32 |x - loc|, not eps32 r. For a pair of adjacent
# cells of width C, |x_i - loc_j| < 2 C and |x_j - loc_j| < C on each axis, so one component of dx
# is off by at most 3 C eps32 and the vector by sqrt(3) times that.
POSITION_ROUNDING_COEFF = 3.0 * math.sqrt(3.0)
KERNEL_NAME = "Wendland C2"
TIME_RTOL = 1e-5
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
        scheme = f["HydroScheme"].attrs if "HydroScheme" in f else {}
        kernel = (
            np.asarray(scheme["Kernel function"]).ravel()[0]
            if "Kernel function" in scheme
            else None
        )
        out["kernel_name"] = kernel.decode() if isinstance(kernel, bytes) else kernel
        cells = f["Cells/Meta-data"].attrs if "Cells/Meta-data" in f else {}
        out["cell_width"] = (
            float(np.max(np.asarray(cells["size"], float))) if "size" in cells else None
        )
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
    params: dict, rep: Report, startup_log: Optional[str]
) -> None:
    """Fail the run as void if the switch-on fixture parameters are not set.

    SWIFT does not write the transported LW shielding option to used_parameters.yml: it
    sets it from ISRF_H2_line_absorption, passes it to the library, and prints the line
    h2_transported_lw_rate_shielding in its start-up log, which this check reads.
    """
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
    if startup_log is None:
        raise GateExit(EXIT_MISSING, "the check needs --startup-log")
    try:
        text = open(startup_log).read()
    except OSError as err:
        raise GateExit(EXIT_MISSING, f"start-up log {startup_log} is missing: {err}")
    if not re.search(r"h2_transported_lw_rate_shielding = 1\b", text):
        raise GateExit(
            EXIT_VOID,
            "void run, the SWIFT start-up line h2_transported_lw_rate_shielding = 1 is absent",
        )
    rep.info(
        "premise parameters",
        "HII, winds, SN off, unlimited speed, SWIFT start-up line reports background_only",
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
    This assumes that every gas particle sits on one time bin: a particle on a shorter bin
    forwards the light more than once per base step, and the gate does not check it.

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
    mean. The kernel mean of LWLineTransmissions is a mean of the neighbours' own
    ratios weighted by their LW energy, so it differs from f(n_H2 r) by the averaging;
    the bars on the own ratio do not carry that offset.
    """
    own = last["LWOwnSpecificEnergies"]
    ratio = last["LWOwnDissociatingSpecificEnergies"] / np.where(own > 0, own, 1.0)
    return own, ratio


def h2_range_premise(
    first: dict,
    last: dict,
    mask: np.ndarray,
    rep: Report,
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
    ranges = h2_region_ranges(first, last, [("", mask)], rep)
    return ranges[0][1], ranges[0][2]


def h2_region_ranges(
    first: dict, last: dict, regions: List[Tuple[str, np.ndarray]], rep: Report
) -> List[Tuple[str, float, float, float, float]]:
    """Apply the n_H2 range premise per region.

    Regions are for domains whose gas differs by design (a slab and clear gas), where
    one range would carry no meaning.

    Returns (label, n_min, n_max, lo_ref, hi_ref) per region, n_min and n_max over the
    reference and the last snapshot, lo_ref and hi_ref over the reference alone.

    Raises
    ------
    GateExit
        EXIT_VOID if a particle of a region left the reference range of that region.
    """
    n_first, n_last = h2_column_density(first), h2_column_density(last)
    out = []
    for label, mask in regions:
        tag = f" ({label})" if label else ""
        if not mask.any():
            raise GateExit(EXIT_VOID, f"void run, the n_H2 domain{tag} holds no gas")
        nf, nl = n_first[mask], n_last[mask]
        lo_ref, hi_ref = float(nf.min()), float(nf.max())
        n_out = int(
            np.sum(
                (nl < lo_ref * (1.0 - FROZEN_H2_REL))
                | (nl > hi_ref * (1.0 + FROZEN_H2_REL))
            )
        )
        if n_out:
            raise GateExit(
                EXIT_VOID,
                f"void run, H2 changed: n_H2 of {n_out} probe-domain particles{tag} left the "
                f"reference range {lo_ref:.5g} to {hi_ref:.5g} cm^-3 by more than "
                f"{FROZEN_H2_REL} (last range {float(nl.min()):.5g} to {float(nl.max()):.5g})",
            )
        n_min, n_max = min(lo_ref, float(nl.min())), max(hi_ref, float(nl.max()))
        rep.info(
            f"n_H2 range premise{tag}",
            f"all domain particles inside the reference range; bar range {n_min:.5g} to {n_max:.5g} cm^-3",
        )
        out.append((label, n_min, n_max, lo_ref, hi_ref))
    return out


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


def kernel_weight_with_error(
    u: np.ndarray, cell_over_h: float
) -> Tuple[np.ndarray, np.ndarray]:
    """Return the code's kernel W(u) for h = 1 and a bound on its float error.

    W is the piecewise polynomial of the code evaluated in double. The error bound is
    the sum of three float effects: the Horner rounding, 10 eps32 sum_k |c_k| |t|^(5-k)
    (the standard bound 2 n eps for degree n = 5); the rounding of t = u - origin, which
    is exact near the edge, so only |dW/dt| eps32 |t|; and the rounding of u itself,
    |dW/du| eps32 (POSITION_ROUNDING_COEFF C / h + ARGUMENT_ROUNDING_EPS u), with C the
    cell width. Their sum is large relative to W only at the edge of the support, where W
    vanishes as (gamma - u)^4.

    Returns
    -------
    (np.ndarray, np.ndarray)
        W and the absolute error bound, both zero at and beyond the code's support.
    """
    u = np.asarray(u, float)
    index = np.minimum(
        (np.minimum(u, KERNEL_CODE_GAMMA) * (4.0 / KERNEL_GAMMA)).astype(int), 4
    )
    index = np.where(u >= KERNEL_CODE_GAMMA, 4, index)
    c = KERNEL_POLY_COEFFS[index]
    t = u - KERNEL_POLY_ORIGIN[index]
    w = c[:, 0].copy()
    dw = np.zeros_like(u)
    mag = np.abs(c[:, 0])
    for k in range(1, 6):
        dw = dw * t + w
        w = w * t + c[:, k]
        mag = mag * np.abs(t) + np.abs(c[:, k])
    w = np.maximum(w, 0.0)
    err = FLOAT32_EPS * (
        10.0 * mag
        + np.abs(dw) * np.abs(t)
        + np.abs(dw)
        * (POSITION_ROUNDING_COEFF * cell_over_h + ARGUMENT_ROUNDING_EPS * u)
        + w
    )
    inside = u < KERNEL_CODE_GAMMA
    return np.where(inside, w, 0.0), np.where(inside, err, 0.0)


def kernel_neighbour_weights(
    snap: dict,
) -> Iterator[Tuple[int, np.ndarray, np.ndarray, np.ndarray]]:
    """Yield (i, j, w, rel) for every particle i.

    j are the neighbours, w = m_j W(r_ij / h_i) the weights of the force loop and rel the
    relative error bound of each weight (see kernel_weight_with_error; the cell width is
    snap["cell_width"]). The list includes
    i itself, which carries the self term m_i W(0) that the code adds when it closes the
    mean. A weight of zero (the edge of the support) carries a relative error of zero; a
    neighbour within a few eps32 of the edge that the code counts and this list does not, or the
    other way round, has a weight below 1e-26 of W(0) and is not accounted for.
    """
    from scipy.spatial import cKDTree

    box = snap["box"]
    cell = snap["cell_width"]
    pos = np.mod(snap["Coordinates"], box)
    h, m = snap["SmoothingLengths"], snap["Masses"]
    tree = cKDTree(pos, boxsize=box)
    ngb = tree.query_ball_point(pos, KERNEL_GAMMA * h * (1.0 + 1e-6))
    for i, found in enumerate(ngb):
        j = np.asarray(found)
        d = pos[j] - pos[i]
        d -= box * np.round(d / box)
        w, err = kernel_weight_with_error(
            np.sqrt((d * d).sum(axis=1)) / h[i], cell / h[i]
        )
        rel = np.where(w > 0.0, err / np.where(w > 0.0, w, 1.0), 0.0)
        yield i, j, m[j] * w, rel


def step_premises(path: Optional[str], snap: dict, rep: Report) -> float:
    """Return the length of the last step; void unless every particle was updated in it.

    The kernel mean is closed in the last step, so it equals the mean of the own fields
    in the snapshot only if all particles were active in that step and the positions
    moved by less than a float rounding of the kernel since then.

    Raises
    ------
    GateExit
        With EXIT_VOID if the file is not given, no step ends at the snapshot time, a
        step did not update every particle or the positions moved; with EXIT_MISSING if
        the file is unreadable.
    """
    if not path:
        raise GateExit(
            EXIT_VOID,
            "the kernel-mean check needs --timesteps (timesteps.txt of the run)",
        )
    rows = []
    try:
        with open(path) as f:
            for line in f:
                tok = line.split()
                if tok and not tok[0].startswith("#") and len(tok) >= 8:
                    rows.append((float(tok[1]), float(tok[4]), int(tok[7])))
    except (OSError, ValueError) as err:
        raise GateExit(EXIT_MISSING, f"timesteps file {path} is unreadable: {err}")
    t_snap = snap["time"]
    hit = [
        r for r in rows if abs(r[0] - t_snap) <= TIME_RTOL * max(abs(t_snap), 1e-300)
    ]
    if not hit:
        raise GateExit(
            EXIT_VOID, f"no step of {path} ends at the snapshot time {t_snap:.7g}"
        )
    _, dt, updates = hit[-1]
    n_gas = len(snap["ParticleIDs"])
    if updates != n_gas:
        raise GateExit(
            EXIT_VOID,
            f"the step that ends at the snapshot updated {updates} of {n_gas} gas particles",
        )
    drift = float(np.max(np.linalg.norm(snap["Velocities"], axis=1))) * dt
    h_min = float(np.min(snap["SmoothingLengths"]))
    if drift > FLOAT32_EPS * h_min:
        raise GateExit(
            EXIT_VOID,
            f"positions moved up to {drift / h_min:.3g} h in the last step, above one float rounding",
        )
    rep.info(
        "kernel-mean premises",
        f"all {n_gas} particles updated in the last step; largest displacement {drift / h_min:.3g} h",
    )
    return dt


def check_kernel_mean(
    last: dict,
    rep: Report,
    probes: np.ndarray,
    exp: np.ndarray,
    bar_own: np.ndarray,
) -> None:
    """Check the snapshot kernel means of the LW fields in three tiers.

    The tiers read the own fields of the snapshot, so the own LW energy E_j itself is
    not tested here. The snapshot must hold the Wendland C2 kernel (else the run is void)
    and the cell width of Cells/Meta-data/size, which bounds the float rounding of the
    positions in the force loops.

    The float LW sum of the code flushes a product w E_j below the smallest normal float
    (FLT_MIN) to zero, or keeps it, depending on the order of the sum and on fused
    multiply-add; the dissociating sum, in double, loses nothing. A product within its own
    relative error of FLT_MIN is undecided. The LW sum lies between the sum of the products
    certainly kept and the sum of all products, and the offline value is measured from the
    window that this range opens. The edge terms are bounded over every set of dropped
    products. A particle whose terms are all below FLT_MIN must hold no LW and no
    dissociating energy.

    Tier 1 (identity): T_offline is the code's kernel mean, sum_j m_j W(r_ij / h_i) D_j
    over sum_j m_j W(r_ij / h_i) E_j with the self term, of the own fields in the
    snapshot (D_j clamped to E_j). It must agree with the snapshot to (n_ngb + 2) eps32.
    A particle with no lit neighbour must have no LW energy. Each bar adds the kernel-edge
    term sum_j max |b_j - a_j| rel_j, with b_j the share of
    neighbour j in the numerator, a_j its share of the weight sum (the largest value over the
    possible sums of the terms kept) and rel_j the relative
    float error of its weight (kernel_weight_with_error): the weights of the neighbours near
    the edge of the support carry a large relative error, and one lit neighbour there can set
    the mean. The weight sum of the mean cancels in T and is tested by tier 3.

    Tier 3 (means): E_mean and D_mean, which the cooling reads, against the offline sum
    w E / sum w and min(sum w D, sum w E) / sum w with the self term. The bar is the float
    sums, (2 n_ngb + 4) eps32 for E_mean (numerator, weight sum, division) and the same
    plus (n_ngb + 3) eps32 for D_mean, plus the kernel-edge term. A change of the weight
    sum, the weights or the fields of 0.5 per cent or more of E_mean or D_mean fails
    wherever the bar is below that. The run is void if no particle has such a bar.

    Tier 2 (physics): T_snapshot against K[f], the same weighted mean of the prediction
    f(n_H2 r) of each neighbour. T / K is a convex combination of the ratios t_j / f_j,
    so the bar is the largest own bar over the lit neighbours plus (n_ngb + 2) eps32.
    It skips the particles that lose more than the float term to the flush, because the
    convex combination does not hold for them.

    Raises
    ------
    GateExit
        With EXIT_MISSING if the snapshot lacks the kernel name or the cell width; with
        EXIT_VOID if the kernel is not Wendland C2 or no particle has a tier-3 bar below
        0.5 per cent.
    """
    if last.get("kernel_name") is None or last.get("cell_width") is None:
        raise GateExit(
            EXIT_MISSING,
            "the snapshot lacks HydroScheme/Kernel function or Cells/Meta-data/size, which the "
            "kernel-mean check reads",
        )
    if last["kernel_name"] != KERNEL_NAME:
        raise GateExit(
            EXIT_VOID,
            f"the snapshot kernel is {last['kernel_name']!r}; the kernel-mean check holds the "
            f"polynomial of {KERNEL_NAME!r} only",
        )
    energy = last["LWOwnSpecificEnergies"]
    dissoc = last["LWOwnDissociatingSpecificEnergies"]
    t_snap = last["LWLineTransmissions"]
    n = len(energy)
    t_win = np.full((2, n), np.nan)
    e_win = np.full((2, n), np.nan)
    d_win = np.full((2, n), np.nan)
    kmean = np.full(n, np.nan)
    bar2 = np.full(n, np.nan)
    n_ngb = np.zeros(n)
    flushed = np.zeros(n)
    gone = np.zeros(n, dtype=bool)
    unlit = np.zeros(n, dtype=bool)
    edge_t = np.zeros(n)
    edge_e = np.zeros(n)
    edge_d = np.zeros(n)
    for i, j, w, rel in kernel_neighbour_weights(last):
        n_ngb[i] = len(j)
        terms = w * energy[j]
        total = terms.sum()
        if total <= 0.0:
            unlit[i] = True
            continue
        # A product below FLT_MIN can be flushed to zero in the float LW sum, depending on the
        # order of the sum and on fused multiply-add. A product within its own relative error of
        # FLT_MIN is undecided. The LW sum lies between the sum of the products certainly kept
        # and the sum of all products; the dissociating sum, in double, loses nothing.
        slack = rel + 2.0 * FLOAT32_EPS
        kept = terms[terms * (1.0 - slack) >= FLOAT32_TINY].sum()
        flushed[i] = terms[terms < FLOAT32_TINY].sum() / total
        gone[i] = not np.any(terms * (1.0 + slack) >= FLOAT32_TINY)
        weight_sum = w.sum()
        d_terms = w * dissoc[j]
        d_total = d_terms.sum()
        e_win[:, i] = np.array([kept, total]) / weight_sum
        d_win[:, i] = np.minimum(d_total, [kept, total]) / weight_sum
        t_win[:, i] = np.minimum(
            1.0, d_total / np.array([total, kept if kept > 0.0 else np.inf])
        )
        share_w = w / weight_sum
        share_d = d_terms / d_total if d_total > 0.0 else np.zeros_like(w)
        # The share of term j in the LW sum is terms_j s with s between 1 / total and 1 / kept, or
        # zero if the term is dropped. |b - share| is convex in s, so its maximum over the
        # possible sums is at an end of the range, or at the dropped value |b|.
        sure = terms * (1.0 - slack) >= FLOAT32_TINY
        share_lo = terms / total
        share_hi = terms / kept if kept > 0.0 else np.ones_like(w)

        def worst_gap(b):
            gap = np.maximum(np.abs(b - share_lo), np.abs(b - share_hi))
            return np.where(sure, gap, np.maximum(gap, np.abs(b)))

        edge_t[i] = float(np.sum(worst_gap(share_d) * rel))
        edge_e[i] = float(np.sum(worst_gap(share_w) * rel))
        edge_d[i] = float(np.sum(np.abs(share_d - share_w) * rel))
        if probes[i]:
            kmean[i] = float(terms @ exp[j]) / total
            bar2[i] = bar_own[j][energy[j] > 0.0].max()
    float_term = (n_ngb + 2.0) * FLOAT32_EPS
    bar_t = float_term + edge_t

    e_snap = last["LWSpecificEnergies"]
    d_snap = last["LWDissociatingSpecificEnergies"]
    stray = unlit & (e_snap != 0.0)
    gone_bad = gone & ((e_snap != 0.0) | (d_snap != 0.0))
    flushable = ~unlit & (flushed > float_term)
    compared = ~unlit & ~gone
    if not compared.any():
        raise GateExit(
            EXIT_VOID, "no particle has a lit neighbour that survives the flush"
        )
    dark = e_snap == 0.0

    def window_deviation(value, window):
        """Return the log distance of value from [window[0], window[1]], zero inside it."""
        with np.errstate(divide="ignore", invalid="ignore"):
            below = np.where(
                value < window[0],
                np.log(window[0] / np.where(value > 0.0, value, 1.0)),
                0.0,
            )
            below = np.where((value <= 0.0) & (window[0] > 0.0), np.inf, below)
            above = np.where(value > window[1], np.log(value / window[1]), 0.0)
        dev = np.maximum(below, above)
        return np.where(np.isfinite(value), dev, np.inf)

    dev1 = np.where(dark & (e_win[0] == 0.0), 0.0, window_deviation(t_snap, t_win))
    over1 = int(np.sum(compared & (dev1 > bar_t)))
    margin1 = np.where(compared, dev1 / bar_t, -np.inf)
    worst1 = int(np.argmax(margin1))
    rep.check(
        "kernel mean against the kernel mean of the own fields (tier 1)",
        over1 == 0 and not stray.any() and not gone_bad.any(),
        f"{over1} of {int(compared.sum())} particles outside |ln(T/T_offline)| <= (n_ngb + 2) eps32 plus the kernel-edge term, "
        f"measured from the window that terms below FLT_MIN open (any change of T; "
        f"the weight sum of the mean is tested by tier 3); "
        f"worst margin {float(margin1[worst1]):.3g} of its bar (deviation {float(dev1[worst1]):.3g}, bar {float(bar_t[worst1]):.3g}), "
        f"median bar {float(np.median(bar_t[compared])):.3g}; "
        f"{int(stray.sum())} of {int(unlit.sum())} particles without a lit neighbour hold LW energy; "
        f"{int(np.sum(flushable & compared))} particles can lose more than the float term to terms below FLT_MIN, "
        f"{int(gone.sum())} lose all their terms and {int(gone_bad.sum())} of those hold LW energy",
    )

    bar_e = float_term + (n_ngb + 2.0) * FLOAT32_EPS + edge_e
    bar_d = bar_e + (n_ngb + 3.0) * FLOAT32_EPS + edge_d
    dev_e = window_deviation(e_snap, e_win)
    dev_d = window_deviation(d_snap, d_win)
    over_e = compared & (dev_e > bar_e)
    over_d = compared & (dev_d > bar_d)
    margin_e = np.where(compared, dev_e / bar_e, -np.inf)
    margin_d = np.where(compared, dev_d / bar_d, -np.inf)
    worst_e, worst_d = int(np.argmax(margin_e)), int(np.argmax(margin_d))
    resolve = int(np.sum(compared & (bar_e < 5e-3)))
    if resolve == 0:
        raise GateExit(
            EXIT_VOID,
            "no particle has a tier-3 bar below 0.5 per cent, so a change of the weight sum, the weights "
            "or the fields of that size cannot be seen",
        )
    rep.check(
        "kernel means of the LW and dissociating energies against their offline means (tier 3)",
        not over_e.any() and not over_d.any(),
        f"{int(over_e.sum())} of {int(compared.sum())} particles outside |ln(E_mean/E_offline)| and "
        f"{int(over_d.sum())} outside |ln(D_mean/D_offline)|; offline mean = sum w E / sum w with the self term, "
        f"measured from the window that terms below FLT_MIN open, "
        f"bar = float sums + sum_j |b_j - a_j| (error of w_j) with the kernel-edge error of each weight; "
        f"worst margin E {float(margin_e[worst_e]):.3g} of its bar (deviation {float(dev_e[worst_e]):.3g}, bar {float(bar_e[worst_e]):.3g}), "
        f"worst margin D {float(margin_d[worst_d]):.3g} (deviation {float(dev_d[worst_d]):.3g}, bar {float(bar_d[worst_d]):.3g}); "
        f"{resolve} of {int(compared.sum())} particles have a bar below 0.5 per cent (a change of 0.5 per cent or more of "
        f"E_mean or D_mean fails there, whatever its cause: the weight sum, the weights or the fields)",
    )

    use = probes & ~flushable & np.isfinite(kmean) & compared
    bar = bar2 + float_term
    dev2 = np.abs(np.log(np.maximum(t_snap, 1e-300) / kmean))
    over2 = use & (dev2 > bar)
    rep.check(
        "kernel mean against the kernel mean of f(n_H2 r) (tier 2)",
        not over2.any(),
        f"{int(over2.sum())} of {int(use.sum())} probes outside |ln(T/K[f])| <= largest own bar of the lit "
        f"neighbours + (n_ngb + 2) eps32; worst margin {float(np.max((dev2 / bar)[use])):.3g} of the bar; "
        f"median bar {float(np.median(bar[use])):.3g}",
    )


def gate_starlit(args: argparse.Namespace) -> int:
    rep = Report()
    fix = load_json(args.fixture, "geometry file")
    first = load_snapshot(args.first, REQUIRED_GAS)
    fields = (
        REQUIRED_GAS + REQUIRED_LW + REQUIRED_OWN + REQUIRED_COLUMN + REQUIRED_KERNEL
    )
    last = load_snapshot(args.last, fields)
    require_finite("last snapshot", **{k: last[k] for k in fields})
    check_box_faces(first)
    params = read_params(args.params)
    calibrate_temperatures(first, last, params, rep)
    check_premise_parameters(params, rep, args.startup_log)
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
    # The carried column is a sum of the receivers' n_H2 times paths that add up to r.
    col_field = last[PARITY_COLUMN] / last["unit_length"] ** 2  # stored in U_L^-2
    col_lo, col_hi = n_range[0] * r * (1.0 - float_term), n_range[1] * r * (
        1.0 + float_term
    )
    col_out = lit & ((col_field < col_lo) | (col_field > col_hi))
    rep.check(
        "carried column field between n_min r and n_max r",
        not col_out.any(),
        f"{int(col_out.sum())} of {int(lit.sum())} lit probes outside "
        f"[n_min r (1 - e), n_max r (1 + e)], e the float term; "
        f"ln(N_field / (n_H2 r)) median {float(np.median(np.log(np.maximum(col_field[lit], 1e-300) / column[lit]))):.4g}",
    )
    lw = last["LWSpecificEnergies"]
    mean_ratio = last["LWDissociatingSpecificEnergies"] / np.where(lw > 0, lw, 1.0)
    lit_mean = read & (lw > 0.0)
    if lit_mean.any():
        rep.reading(
            "kernel-mean ratio LWLineTransmissions",
            f"median ln(ratio/f) {float(np.median(np.log(np.maximum(mean_ratio[lit_mean], 1e-300) / exp[lit_mean]))):.4g} "
            f"over {int(lit_mean.sum())} probes; it exceeds f where brighter neighbours are "
            f"nearer the star, a property of the averaging (checked against the kernel mean of f below)",
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

    step_premises(args.timesteps, last, rep)
    probes = read & (last["LWSpecificEnergies"] > 0.0)
    check_kernel_mean(last, rep, probes, exp, bar)
    return rep.exit_code()


LOG_IDENTITY = (
    ("code version", r"^\s*Version\s*:.*$"),
    ("code revision", r"^\s*Revision:.*$"),
    ("compiler", r"^\s*Compiler:.*$"),
    ("geajay library", r"geajay version = .*$"),
    ("slot frame", r"octant frame rotation \(row-major\)\s*=.*$"),
)
SLOTS_HANDLED = 12


def read_startup_log(path: str) -> Dict[str, str]:
    """Return the lines of a swift log that identify the build and the slot frame.

    Raises
    ------
    GateExit
        EXIT_MISSING if the log or one of the lines is absent.
    """
    try:
        text = open(path).read()
    except OSError as err:
        raise GateExit(EXIT_MISSING, f"start-up log {path} is missing: {err}")
    out = {}
    for name, pattern in LOG_IDENTITY:
        found = re.search(pattern, text, re.M)
        if found is None:
            raise GateExit(EXIT_MISSING, f"start-up log {path} has no {name} line")
        out[name] = found.group(0).strip()
    return out


def slot_frame(line: str) -> np.ndarray:
    """Return the 3x3 slot frame (row-major) from its start-up log line."""
    vals = [float(v) for v in line.split("=", 1)[1].replace(",", " ").split()]
    if len(vals) != 9:
        raise GateExit(
            EXIT_MISSING, f"slot frame line holds {len(vals)} numbers, not 9"
        )
    return np.array(vals).reshape(3, 3)


def packet_slot(frame: np.ndarray, direction: np.ndarray, n_slots: int) -> np.ndarray:
    """Return the reception slot of each direction, as radiation_packet_slot for 12 slots."""
    if n_slots != SLOTS_HANDLED:
        raise GateExit(
            EXIT_MISSING,
            f"the check handles {SLOTS_HANDLED} reception slots, the build has {n_slots}",
        )
    q = np.atleast_2d(direction) @ frame.T
    m = np.argmin(np.abs(q), axis=1)
    i = np.where(m == 0, 1, 0)
    j = np.where(m == 2, 1, 2)
    r = np.arange(len(q))
    return 4 * m + 2 * (q[r, i] >= 0.0) + (q[r, j] >= 0.0)


def image_displacements(pos: np.ndarray, star: np.ndarray, box: float):
    """Return the vectors and lengths from the 27 nearest periodic images of a star to pos."""
    shifts = box * np.array(list(itertools.product((-1.0, 0.0, 1.0), repeat=3)))
    d = pos[:, None, :] - (star[None, None, :] + shifts[None, :, :])
    return d, np.linalg.norm(d, axis=2)


def star_particles(path: str) -> Tuple[np.ndarray, np.ndarray]:
    """Return the coordinates and the IDs of the star particles of a snapshot."""
    try:
        with h5py.File(path, "r") as f:
            if "PartType4" not in f:
                raise GateExit(EXIT_MISSING, f"{path} has no PartType4")
            return (
                np.asarray(f["PartType4/Coordinates"][:], float).reshape(-1, 3),
                np.asarray(f["PartType4/ParticleIDs"][:], dtype=np.int64),
            )
    except OSError as err:
        raise GateExit(EXIT_MISSING, f"snapshot {path} cannot be opened: {err}")


def previous_snapshot(path: str) -> str:
    """Return the path of the snapshot written before the given one."""
    found = re.search(r"snapshot_(\d+)\.hdf5$", path)
    if found is None or int(found.group(1)) < 1:
        raise GateExit(
            EXIT_MISSING, f"{path} is not a numbered snapshot after the first"
        )
    num = found.group(1)
    return path[: found.start(1)] + f"{int(num) - 1:0{len(num)}d}.hdf5"


def run_last_snapshot(run_dir: str, like: str) -> str:
    """Return the snapshot of a run directory with the same file name as the reference."""
    path = os.path.join(run_dir, "snap", os.path.basename(like))
    if not os.path.exists(path):
        raise GateExit(
            EXIT_VOID,
            f"void run, {run_dir} has no {os.path.basename(like)} like the two-star run",
        )
    return path


def align_to_ids(snap: dict, ids: np.ndarray, label: str) -> None:
    """Reorder the per-particle arrays of a snapshot to the order of ids."""
    order = np.argsort(snap["ParticleIDs"])
    where = np.minimum(np.searchsorted(snap["ParticleIDs"][order], ids), len(order) - 1)
    if len(snap["ParticleIDs"]) != len(ids) or not np.array_equal(
        snap["ParticleIDs"][order][where], ids
    ):
        raise GateExit(EXIT_VOID, f"void run, {label} holds other particles")
    idx = order[where]
    for key, val in snap.items():
        if isinstance(val, np.ndarray) and val.ndim >= 1 and len(val) == len(ids):
            snap[key] = val[idx]


def ray_column(
    star: np.ndarray,
    pos: np.ndarray,
    box: float,
    slab: Tuple[float, float],
    n_clear: float,
    n_slab: float,
    reach: float,
) -> Tuple[np.ndarray, np.ndarray]:
    """Return the straight-ray column from a star (nearest image) and the face variation.

    The gas is n_clear outside the slab (periodic in x) and n_slab inside it. Positions are
    in the length unit of the box. The column is in n times that unit. The variation is
    the jump |n_slab - n_clear| times the number of slab faces within reach of the ray in
    x: the SPH density spreads a face over the kernel support, so a face closer than
    reach to the ray is crossed by the smoothed field as well. The size of one jump per
    face follows from the smoothing reach; it is not a proven bound for a slanted ray.
    The reference does not call the shielding function.
    """
    d = pos - star[None, :]
    d -= box * np.round(d / box)
    length = np.linalg.norm(d, axis=1)
    x_end = star[0] + d[:, 0]
    lo, hi = np.minimum(star[0], x_end), np.maximum(star[0], x_end)
    overlap, faces = np.zeros(len(pos)), np.zeros(len(pos))
    star_in_slab = False
    for m in range(-2, 3):
        a, b = slab[0] + m * box, slab[1] + m * box
        overlap += np.clip(np.minimum(hi, b) - np.maximum(lo, a), 0.0, None)
        faces += ((a > lo - reach) & (a < hi + reach)).astype(float)
        faces += ((b > lo - reach) & (b < hi + reach)).astype(float)
        star_in_slab |= a <= star[0] <= b
    dx = np.abs(d[:, 0])
    fraction = np.where(
        dx > 0.0, overlap / np.where(dx > 0.0, dx, 1.0), float(star_in_slab)
    )
    in_slab = length * fraction
    return (
        n_clear * length + (n_slab - n_clear) * in_slab,
        abs(n_slab - n_clear) * faces,
    )


def round_off_bar(steps: float, k_terms: int, p_power: float, delta_h: float) -> float:
    """Return n (k 2^-24 + p delta_h), the round-off bar of a field that is a sum of shares."""
    return steps * (k_terms * FLOAT32_EPS + p_power * delta_h)


def smoothing_gap(h_x: np.ndarray, h_y: np.ndarray) -> float:
    """Return the largest relative smoothing-length gap between two runs."""
    return float(np.max(np.abs(h_x - h_y) / np.maximum(h_x, h_y)))


def round_off_gaps(
    x: np.ndarray, y: np.ndarray, held: np.ndarray, ngb: List[List[int]], bar: float
) -> Dict[str, float]:
    """Return how far two fields differ in the round-off sense of the parity check.

    Particles in held are held to bar times their own value. The others (light front and
    dark) are held to bar times the largest value in their kernel.
    """
    scale = np.maximum(np.abs(x), np.abs(y))
    gap = np.abs(x - y)
    own = np.where(scale > 0, gap / np.where(scale > 0, scale, 1), 0.0)
    local = np.array([scale[n].max() for n in ngb])
    rest = ~held
    worst = np.where(local > 0, gap / np.where(local > 0, local, 1), 0.0)
    return {
        "n_held": int(held.sum()),
        "over_held": int(np.sum(held & (gap > bar * scale))),
        "held_own": float(own[held].max(initial=0.0)),
        "n_rest": int(rest.sum()),
        "over_rest": int(np.sum(rest & (gap > bar * local))),
        "rest_local": float(worst[rest].max(initial=0.0)),
        "rest_own": float(own[rest].max(initial=0.0)),
    }


def gate_twosource(args: argparse.Namespace) -> int:
    """Check two stars that never share a reception slot, with one-star runs as references.

    Four runs share one binary: star A alone, star B alone, star A alone again, and both.
    Where the light of the two stars never sits in one slot of a receiver, the field of
    the two-star run is the sum of the one-star fields. Four checks follow.

    1. Premises, read from the runs: same build, parameters and slot frame; the stars
       where the fixture says; the light of the two stars shares no slot at any particle
       that holds both (all periodic images within the reach of the light, which is the
       largest lit distance of the one-star run plus one hop); two identical runs agree
       within the bar of check 2, which controls the bar.
    2. The own LW energy and its dissociating part of the two-star run equal the sums of
       the one-star fields within the round-off bar of the two-rank check,
       n (k 2^-24 + p delta_h), for gas that held light at the previous snapshot, and
       within the same bar times the largest value in the kernel for the light front
       and dark gas. A field is a sum of positive shares, so reordering the sum moves it
       by at most that. A coupling of the two lights moves it more.
    3. The column of each one-star run lies within the interface bar of the design,
       l_hop (n_H2 + TV), of the straight-ray column through the slab and the clear gas
       (nearest periodic image). This sees an error common to the one-star and the
       two-star runs, which check 2 cancels.
    4. Where both stars are lit, the column of the two-star run is the mean of the one-star
       columns weighted by their dissociating energies, N_AB = (D_A N_A + D_B N_B) / (D_A + D_B),
       because the code stores sum(D_s N_s) / sum(D_s) over the slots. Numerator and
       denominator are sums of positive shares. The bar covers the round-off of the energy
       weights, and the ratio lies within twice that bar. The gap of the slot column
       between the runs is not bounded by check 2; it stayed below 0.02 of the bar on the
       default run. The identity needs ISRF_packet_read_mean_time_Myr = 0, because the
       dissociating energy is a running mean and the column is weighted by the
       instantaneous energies.

    Raises
    ------
    GateExit
        EXIT_VOID if a premise is not met.
    """
    rep = Report()
    fix = load_json(args.fixture, "geometry file")
    for key in (
        "box_kpc",
        "star_a_kpc",
        "star_b_kpc",
        "slab_x_kpc",
        "n_h2_clear_cm3",
        "n_h2_slab_cm3",
    ):
        if key not in fix:
            raise GateExit(EXIT_MISSING, f"geometry file lacks {key}")
    fields = REQUIRED_GAS + REQUIRED_OWN + REQUIRED_COLUMN
    first = load_snapshot(args.first, REQUIRED_GAS)
    ab = load_snapshot(args.last, fields)
    require_finite("two-star snapshot", **{k: ab[k] for k in fields})
    check_box_faces(first)
    ids = ab["ParticleIDs"]
    align_to_ids(first, ids, "snapshot 1")
    runs = {"two-star": ab}
    paths = {"two-star": args.last}
    run_dirs = {
        "star A": args.only_a,
        "star B": args.only_b,
        "star A again": args.repeat_a,
    }
    for label, run_dir in run_dirs.items():
        paths[label] = run_last_snapshot(run_dir, args.last)
        runs[label] = load_snapshot(paths[label], fields)
        align_to_ids(runs[label], ids, label)
        require_finite(
            label, **{k: runs[label][k] for k in REQUIRED_OWN + REQUIRED_COLUMN}
        )
        if runs[label]["time"] != ab["time"]:
            raise GateExit(
                EXIT_VOID,
                f"void run, {label} is at time {runs[label]['time']}, the two-star run at {ab['time']}",
            )
    previous = {}
    for label in ("two-star", "star A"):
        snap = load_snapshot(
            previous_snapshot(paths[label]), ("ParticleIDs", "LWOwnSpecificEnergies")
        )
        align_to_ids(snap, ids, f"{label} previous snapshot")
        previous[label] = snap["LWOwnSpecificEnergies"] > 0.0

    params = read_params(args.params)
    check_premise_parameters(params, rep, args.startup_log)
    if not param_equals(
        params, "GEARFeedback:ISRF_packet_read_mean_time_Myr", "0", "0"
    ):
        raise GateExit(
            EXIT_VOID,
            "void run, ISRF_packet_read_mean_time_Myr is not 0: the column weights are instantaneous",
        )
    logs = {"two-star": read_startup_log(args.startup_log)}
    for label, run_dir in run_dirs.items():
        other = read_params(os.path.join(run_dir, "used_parameters.yml"))
        if other != params:
            diff = sorted(
                k for k in set(other) | set(params) if other.get(k) != params.get(k)
            )
            raise GateExit(
                EXIT_VOID, f"void run, {label} used other parameters: {diff[:5]}"
            )
        logs[label] = read_startup_log(os.path.join(run_dir, "output.log"))
    for name, _ in LOG_IDENTITY:
        if len({log[name] for log in logs.values()}) != 1:
            raise GateExit(
                EXIT_VOID, f"void run, the four runs differ in the {name} line"
            )
    rep.info(
        "same build and parameters",
        "four runs; "
        + "; ".join(f"{v}" for k, v in logs["two-star"].items() if k != "slot frame"),
    )
    frame = slot_frame(logs["two-star"]["slot frame"])
    box = ab["box"]
    star_a = np.asarray(fix["star_a_kpc"], float)
    star_b = np.asarray(fix["star_b_kpc"], float)
    expect = {
        "two-star": [star_a, star_b],
        "star A": [star_a],
        "star A again": [star_a],
        "star B": [star_b],
    }
    star_id = {}
    for label, want in expect.items():
        coords, sids = star_particles(paths[label])
        dist = [np.linalg.norm(coords - w[None, :], axis=1) for w in want]
        if len(coords) != len(want) or any(d.min() >= 1e-5 * box for d in dist):
            raise GateExit(
                EXIT_VOID,
                f"void run, the stars of the {label} run are not at the fixture positions",
            )
        star_id[label] = [int(sids[int(np.argmin(d))]) for d in dist]
    for label, want_ids in (
        ("star A", star_id["two-star"][:1]),
        ("star A again", star_id["two-star"][:1]),
        ("star B", star_id["two-star"][1:]),
    ):
        if star_id[label] != want_ids:
            raise GateExit(
                EXIT_VOID,
                f"void run, the star of the {label} run has ID {star_id[label]}, in the two-star run {want_ids}: "
                "the ID seeds the random frame of the emission pixels",
            )
    rep.info(
        "stars",
        f"the two-star run holds both stars, the one-star runs star A, star B and star A, at the "
        f"fixture positions and with the IDs of the two-star run ({star_id['two-star']})",
    )

    slab = tuple(fix["slab_x_kpc"])
    in_slab = (ab["Coordinates"][:, 0] >= slab[0]) & (
        ab["Coordinates"][:, 0] <= slab[1]
    )
    h2_region_ranges(first, ab, [("slab gas", in_slab), ("clear gas", ~in_slab)], rep)

    own = {k: runs[k]["LWOwnSpecificEnergies"] for k in runs}
    dis = {k: runs[k]["LWOwnDissociatingSpecificEnergies"] for k in runs}
    lit = {k: own[k] > 0.0 for k in runs}
    both = lit["star A"] & lit["star B"] & lit["two-star"]
    if not both.any():
        raise GateExit(EXIT_VOID, "void run, no particle is lit by both stars")

    control = packet_slot(
        frame,
        np.array([[1.0, 0, 0], [2.0, 0, 0], [1.0, 0.01, 0], [-1.0, 0, 0]]),
        args.slots,
    )
    if not (control[0] == control[1] == control[2] != control[3]):
        raise GateExit(
            EXIT_VOID, f"void check, the slot rule fails its control {control.tolist()}"
        )
    pos = ab["Coordinates"]
    h_max = float(ab["SmoothingLengths"].max())
    images, reach = {}, {}
    for label, star in (("star A", star_a), ("star B", star_b)):
        images[label] = image_displacements(pos, star, box)
        dist = images[label][1]
        reach[label] = float(dist.min(axis=1)[lit[label]].max()) + KERNEL_GAMMA * h_max
        if reach[label] >= box:
            raise GateExit(
                EXIT_VOID,
                f"void check, the light of {label} reaches {reach[label]:.4g}, a box is {box:.4g}",
            )
        many = int(np.sum(np.sum(dist <= reach[label], axis=1)[lit[label]] > 1))
        if many:
            raise GateExit(
                EXIT_VOID,
                f"void run, {many} lit particles are in reach of two images of {label}",
            )
    slots = {
        k: np.stack(
            [packet_slot(frame, images[k][0][:, i, :], args.slots) for i in range(27)],
            axis=1,
        )
        for k in images
    }
    near = {k: images[k][1] <= reach[k] for k in images}
    shared = np.zeros(len(pos), bool)
    for i in range(27):
        for j in range(27):
            shared |= (
                near["star A"][:, i]
                & near["star B"][:, j]
                & (slots["star A"][:, i] == slots["star B"][:, j])
            )
    n_shared = int(np.sum(shared & both))
    h_min = float(ab["SmoothingLengths"].min())
    rep.info(
        "slot premise",
        f"reach of the light {reach['star A'] / h_min:.4g} and {reach['star B'] / h_min:.4g} h_min "
        f"(largest lit distance plus gamma h_max); {int(both.sum())} particles lit by both stars; "
        f"{n_shared} hold both lights in one slot",
    )
    if n_shared:
        raise GateExit(
            EXIT_VOID,
            f"void run, {n_shared} particles hold the light of both stars in one slot: "
            "the sum of the one-star runs is not the reference there",
        )

    steps = float(light_hops(first, ab))
    ngb = kernel_neighbours(ab)
    k_terms = max(len(n) for n in ngb)
    p_power = float(params.get("GEARFeedback:ISRF_packet_cone_power", 4.0)) + 4.0
    h = {k: runs[k]["SmoothingLengths"] for k in runs}
    delta_sum = max(
        smoothing_gap(h["two-star"], h["star A"]),
        smoothing_gap(h["two-star"], h["star B"]),
    )
    delta_rep = smoothing_gap(h["star A"], h["star A again"])
    bar_sum = round_off_bar(steps, k_terms, p_power, delta_sum)
    bar_rep = round_off_bar(steps, k_terms, p_power, delta_rep)
    rep.info(
        "round-off bar",
        f"n (k 2^-24 + p delta_h): sum {steps:g} ({k_terms} x 2^-24 + {p_power:g} x {delta_sum:.3g}) "
        f"= {bar_sum:.3g}; repeat {bar_rep:.3g}",
    )
    held_sum, held_rep = previous["two-star"], previous["star A"]
    if not held_sum.any() or not held_rep.any():
        raise GateExit(
            EXIT_VOID, "void run, no particle held light at the previous snapshot"
        )
    for name, key in (("LW energy", own), ("dissociating energy", dis)):
        g = round_off_gaps(key["star A"], key["star A again"], held_rep, ngb, bar_rep)
        rep.info(
            "control of the bar",
            f"two identical star A runs, {name}: largest gap {g['held_own']:.3g} of the own value "
            f"({g['held_own'] / bar_rep:.3g} of the bar), front and dark {g['rest_local']:.3g} of the "
            f"kernel maximum ({g['rest_local'] / bar_rep:.3g} of the bar)",
        )
        if g["over_held"] or g["over_rest"]:
            raise GateExit(
                EXIT_VOID,
                f"void check, two identical runs differ by more than the round-off bar in the {name}: "
                "the bar is unfit",
            )

    for name, key in (("own LW energy", own), ("own dissociating LW energy", dis)):
        g = round_off_gaps(
            key["two-star"], key["star A"] + key["star B"], held_sum, ngb, bar_sum
        )
        rep.check(
            f"{name}: two stars equal the sum of the one-star runs, lit before",
            g["over_held"] == 0,
            f"{g['over_held']} of {g['n_held']} particles over the bar, largest gap over the own value "
            f"{g['held_own']:.3g}, bar {bar_sum:.3g}",
        )
        rep.check(
            f"{name}: two stars equal the sum of the one-star runs, light front and dark",
            g["over_rest"] == 0,
            f"{g['over_rest']} of {g['n_rest']} particles over the bar times the kernel maximum, "
            f"largest {g['rest_local']:.3g}; largest gap over the own value {g['rest_own']:.3g} (read only)",
        )
    d_a, d_b = dis["star A"][both], dis["star B"][both]
    n_a, n_b = runs["star A"][PARITY_COLUMN][both], runs["star B"][PARITY_COLUMN][both]
    n_ab = ab[PARITY_COLUMN][both]
    require_finite("columns of the particles lit by both", a=n_a, b=n_b, two_star=n_ab)
    if not np.all(d_a + d_b > 0.0):
        raise GateExit(
            EXIT_VOID,
            "void run, a particle lit by both has no dissociating energy",
        )
    if not np.all(n_ab > 0.0):
        raise GateExit(
            EXIT_BAR,
            "FAIL: the two-star column is zero at a particle lit by both that has dissociating energy",
        )
    mix = (d_a * n_a + d_b * n_b) / (d_a + d_b)
    gap_col = np.abs(n_ab / mix - 1.0)
    bar_col = 2.0 * bar_sum
    rep.check(
        "column of the two-star run: dissociating-energy-weighted mean of the one-star columns",
        bool(np.all(gap_col <= bar_col)),
        f"{int(np.sum(gap_col > bar_col))} of {int(both.sum())} particles lit by both outside "
        f"|N_AB / ((D_A N_A + D_B N_B) / (D_A + D_B)) - 1| <= 2 x {bar_sum:.3g}; "
        f"largest {float(gap_col.max()):.3g} ({float(gap_col.max()) / bar_col:.3g} of the bar)",
    )
    ratio_ab = dis["two-star"][both] / own["two-star"][both]
    ratio_sum = (dis["star A"][both] + dis["star B"][both]) / (
        own["star A"][both] + own["star B"][both]
    )
    require_finite(
        "own ratios of the particles lit by both",
        two_star=ratio_ab,
        one_star_sum=ratio_sum,
    )
    rep.reading(
        "ln of the two-star ratio over the ratio of the sums",
        f"largest |ln| {float(np.max(np.abs(np.log(ratio_ab / ratio_sum)))):.3g} "
        f"over {int(both.sum())} particles lit by both stars",
    )
    rep.reading(
        "H2 of the two-star run against the one-star runs",
        f"largest |ln| {float(np.max(np.abs(np.log(ab['H2I'] / runs['star A']['H2I'])))):.3g} (star A), "
        f"{float(np.max(np.abs(np.log(ab['H2I'] / runs['star B']['H2I'])))):.3g} (star B); "
        "the dissociation of H2 by the other star is not in the round-off bar",
    )

    cm = ab["unit_length"]
    for label, star in (("star A", star_a), ("star B", star_b)):
        snap = runs[label]
        l_hop = KERNEL_GAMMA * float(snap["SmoothingLengths"].max())
        n_ray, tv = ray_column(
            star,
            snap["Coordinates"],
            box,
            slab,
            fix["n_h2_clear_cm3"],
            fix["n_h2_slab_cm3"],
            l_hop,
        )
        n_ray *= cm
        bar = l_hop * cm * (h2_column_density(snap) + tv)
        n_code = snap[PARITY_COLUMN] / cm**2
        lit_s = lit[label]
        dev = np.abs(n_code - n_ray)
        out = lit_s & (dev > bar)
        seen = int(np.sum(lit_s & (np.abs(1.5 * n_code - n_ray) > bar)))
        rep.check(
            f"column of {label} against the straight ray",
            not out.any(),
            f"{int(out.sum())} of {int(lit_s.sum())} lit particles outside |N - N_ray| <= l_hop (n_H2 + TV); "
            f"worst {float(np.max((dev / bar)[lit_s])):.3g} of the bar; "
            f"ln(N / N_ray) median {float(np.median(np.log(n_code[lit_s] / n_ray[lit_s]))):.3g}",
        )
        rep.info(
            f"column error of 1.5 times, {label}",
            f"outside the bar at {seen} of {int(lit_s.sum())} lit particles",
        )
        if seen == 0:
            raise GateExit(
                EXIT_VOID,
                f"void check, a column of {label} 1.5 times too large is outside the bar nowhere",
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
    check_premise_parameters(params, rep, args.startup_log)
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
    counters = load_json(args.counters, "counters file") if args.counters else {}
    need = (
        ["probes_on_other_rank", "foreign_bookings"]
        if args.mode == "ranks"
        else ["held_lw_slots_at_restart", "refits_after_restart"]
    )
    reduced = args.counters is None
    if not reduced:
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
    if reduced:
        print(
            "REDUCED  no counters file: the run is not shown to exercise the "
            f"comparison ({', '.join(need)} not checked); a pass proves agreement only"
        )
    else:
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
        return EXIT_REDUCED if reduced and not rep.failed else rep.exit_code()

    if "steps" in counters:
        steps = float(counters["steps"])
    else:
        if args.first is None:
            raise GateExit(
                EXIT_MISSING, "ranks mode needs --first or a counters file with steps"
            )
        steps = float(light_hops(load_snapshot(args.first, geom), a))
    ngb = kernel_neighbours(a)
    k_terms = max(len(n) for n in ngb)
    params = read_params(args.params) if args.params else {}
    p_power = float(params.get("GEARFeedback:ISRF_packet_cone_power", 4.0)) + 4.0
    ha, hb = a["SmoothingLengths"][oa], b["SmoothingLengths"][ob]
    delta_h = smoothing_gap(ha, hb)
    bar = round_off_bar(steps, k_terms, p_power, delta_h)
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
        g = round_off_gaps(x, y, held, ngb, bar)
        rep.check(
            f"{key} round-off, lit before",
            g["over_held"] == 0,
            f"{g['over_held']} of {g['n_held']} particles over the bar, "
            f"largest gap over the own value {g['held_own']:.3g}, bar {bar:.3g}",
        )
        rep.check(
            f"{key} round-off, light front and dark",
            g["over_rest"] == 0,
            f"{g['over_rest']} of {g['n_rest']} particles over the bar times "
            f"the kernel maximum, largest {g['rest_local']:.3g}; largest gap over the "
            f"own value {g['rest_own']:.3g} (read only)",
        )
    return EXIT_REDUCED if reduced and not rep.failed else rep.exit_code()


def gate_stale_key(args: argparse.Namespace) -> int:
    """Check that a given GeajayCooling:h2_transported_lw_rate_shielding has no effect.

    The run gave the removed key a value that contradicts ISRF_H2_line_absorption: 1.
    It must start, and the SWIFT start-up line must report background_only, which only the switch sets.
    """
    rep = Report()
    try:
        text = open(args.log).read()
    except OSError as err:
        raise GateExit(EXIT_MISSING, f"log {args.log} is missing: {err}")
    rep.check(
        "run started and stopped cleanly",
        args.exit_status == 0,
        f"exit status {args.exit_status}",
    )
    values = re.findall(r"h2_transported_lw_rate_shielding = (\d+)", text)
    rep.check(
        "SWIFT start-up line reports background_only (1) for the switch on",
        values == ["1"],
        f"start-up values {values}",
    )
    return rep.exit_code()


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    def common(sp: argparse.ArgumentParser, counters: bool = True) -> None:
        sp.add_argument("--first", required=True, help="First snapshot")
        sp.add_argument("--last", required=True, help="Last snapshot")
        sp.add_argument("--fixture", default="fixture.json")
        sp.add_argument("--params", default="used_parameters.yml")
        sp.add_argument(
            "--startup-log",
            default=None,
            help="swift log with the library start-up lines (output.log)",
        )
        if counters:
            sp.add_argument("--counters", default=None)

    for name, fn in (
        ("starlit", gate_starlit),
        ("twosource", gate_twosource),
        ("velocity", gate_velocity),
    ):
        sp = sub.add_parser(name)
        common(sp, counters=name != "twosource")
        if name == "starlit":
            sp.add_argument(
                "--timesteps",
                default=None,
                help="timesteps.txt of the run, for the premises of the kernel-mean check",
            )
        if name == "velocity":
            sp.add_argument(
                "--case", choices=["auto", "reset", "noreset"], default="auto"
            )
        if name == "twosource":
            for flag, what in (
                ("--only-a", "run directory of star A alone"),
                ("--only-b", "run directory of star B alone"),
                ("--repeat-a", "run directory of star A alone, run again"),
            ):
                sp.add_argument(flag, required=True, help=what)
            sp.add_argument(
                "--slots",
                type=int,
                required=True,
                help="RADIATION_PACKET_SLOTS of the build (config.h)",
            )
        sp.set_defaults(fn=fn)
    sp = sub.add_parser("parity")
    sp.add_argument("--a", required=True)
    sp.add_argument("--b", required=True)
    sp.add_argument("--mode", choices=["ranks", "restart"], required=True)
    sp.add_argument(
        "--counters",
        default=None,
        help="counters file of the premises; without it the check is reduced",
    )
    sp.add_argument(
        "--first",
        default=None,
        help="ranks: first snapshot of run a, whose time gives the step count",
    )
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
    sp = sub.add_parser("stale_key")
    sp.add_argument("--log", required=True)
    sp.add_argument("--exit-status", type=int, required=True)
    sp.set_defaults(fn=gate_stale_key)
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
