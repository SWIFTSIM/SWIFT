################################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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
"""Check a star's propagating ISRF in an expanding box against its exact solution.

Run on the ``star_propagation`` configuration: one star of constant band
luminosity ``L_b`` in a periodic, metal-free box (``kappa = 0``), propagation
on at one fixed speed ``c_pin`` for the whole box (``fixed_fraction``), every
particle and the star on the same time bin. Write ``eps = c_pin/c``.

Under these conditions the pairwise transport and dissipation operators cancel
in the box total ``E_b = sum_i m_i u_b,i`` (radiation_propagation_iact.h), so
the totals obey a closed linear system:

    dE_LW/dt = eps L_LW - eps lambda_E(LW) H E_LW
    dE_N/dt  = eps L_LW - eps lambda_N(LW) H E_N
    dE_PE/dt = eps L_PE - eps lambda_E(PE) H E_PE
               + eps [lambda_E(LW) - 1] H E_LW

with ``E_N`` the LW photon-number moment, the lambdas the run's own band-edge
weights (from its log) and ``H`` the Hubble rate. Each parcel of light decays
from its own emission time, so ``E_b`` is a convolution of the source with the
decay. With ``x = ln a(t) - ln a(t')``:

    E_LW(t) = eps L_LW int exp(-eps lambda_E(LW) x) dt'
    E_PE(t) = eps L_PE int exp(-eps lambda_E(PE) x) dt'
              + eps L_LW [lambda_E(LW) - 1]/[lambda_E(LW) - lambda_E(PE)]
                int [exp(-eps lambda_E(PE) x) - exp(-eps lambda_E(LW) x)] dt'

The time integrals run over the light the code has folded into ``u``. The star
injects its step's dose ``L_b Delta_t`` at the end of that step into a
reservoir, which the next gas step folds in (radiation_iact.h and
radiation_isrf.c's radiation_snapshot_part_propagation). A snapshot written
inside step n therefore holds the doses of steps 1 to n-1 in ``u``.

The check compares the measured totals with two models.

- (S1) The discrete recursion of the code's own update, in double:
  ``E_n = exp(-a_n) E_{n-1} + phi(a_n) eps D_{n-1}`` with
  ``a_n = eps lambda H(a_n) dt_n``, ``phi(a) = (1 - exp(-a))/a``, ``H`` at
  the step end (the engine updates the cosmology to the step end before the
  tasks run) and, for PE, the transfer
  ``(lambda_E(LW) - 1)/lambda_E(LW)`` times the LW loss of the step. Its bar
  is float arithmetic only.
- (S2) The continuous solution above. Its bar is the float bar plus the
  discrete-minus-continuous difference, which is computed from the step grid
  and printed. That difference is the time-integration error of the scheme:
  the code holds ``H`` at the step end, while ``int H dt`` over a step is
  exactly ``d ln a``, and folds each dose in one step after its emission. Both
  are first order in the step's ``d ln a``.

Both residuals are relative to ``E_noH``, the folded dose times ``eps``. The
Hubble signal, ``(E_noH - E_cont)/E_noH``, is printed per band. A cosmological
leg FAILS when the LW signal does not exceed its S2 bar, because such a leg
cannot tell an expanding box from a static one. At this fixture's band ratio
the PE signal nearly cancels: its own decay and its gain from LW have opposite
signs. The PE signal is reported, not required.

The float bar is a sum of worst-case bounds, each derived from the code and
none from a measurement:

- the star's step narrowed to a float (radiation_iact.h,
  feedback_star_store_timestep), the dose stored as a float into the
  reservoir, and the reservoir divided by a float step when it is drained:
  one float32 unit roundoff each;
- the weight-sum closure and the reconstruction terms of the injection
  identity, as `injection_identity_bar` in isrf_cosmology_check.py derives
  them, with the star's neighbour count counted inside its kernel radius;
- the half unit in the fifth significant digit with which the log prints each
  lambda, propagated through the Hubble signal;
- the largest relative mismatch between SWIFT's clock (snapshot ``Time``) and
  the Friedmann time the models integrate, since the dose of the code is
  ``L_b`` times SWIFT's own steps.

One float term has no bound computable from a snapshot: the float
accumulation of the pairwise transport and dissipation sums, which cancel in
exact arithmetic. The non-cosmological twin runs the same transport with
``H = 0``, where both models reduce to ``E_noH``; it must pass the same bar.
"""

import argparse
import glob
import os
import sys
from typing import Dict, List, Tuple

import h5py
import numpy as np
import yaml

from isrf_cosmology_check import (
    FLOAT32_EPS,
    KERNEL_GAMMA_DEFAULT,
    load_run,
    read_band_edge_weights,
    read_c_hyp_scheme,
    read_isrf_propagation,
)

U32 = FLOAT32_EPS / 2.0
C_LIGHT_CGS = 2.99792458e10
# Gauss-Legendre nodes per step for the Friedmann and convolution integrals.
GAUSS_NODES, GAUSS_WEIGHTS = np.polynomial.legendre.leggauss(8)
# A snapshot this close to a step boundary, as a fraction of the step, cannot
# be placed in a step from its double-precision time, and is skipped.
STEP_BOUNDARY_GUARD = 1e-6
BANDS = ("PE", "LW", "N")
U_KEYS = {"PE": "u_PE", "LW": "u_LW", "N": "n_LW"}
BAND_NAMES = {"PE": "PE energy", "LW": "LW energy", "N": "LW photon number"}


def parse_options() -> argparse.Namespace:
    """Parse the command line."""
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "-s", "--snapshots", required=True, help="Snapshot glob of the run"
    )
    parser.add_argument(
        "--kernel-gamma",
        type=float,
        default=KERNEL_GAMMA_DEFAULT,
        help="Kernel radius over smoothing length (default: Wendland C2)",
    )
    parser.add_argument(
        "--table", action="store_true", help="Print the residuals at every snapshot"
    )
    parser.add_argument(
        "--no-expansion-in-models",
        action="store_true",
        help="Drop the Hubble terms from both models. A cosmological leg must "
        "then FAIL: this is the check's own wrong-input test",
    )
    return parser.parse_args()


def run_directory(pattern: str) -> str:
    """Return the run directory holding ``snap/``."""
    return os.path.dirname(os.path.dirname(sorted(glob.glob(pattern))[0]))


def read_parameters(pattern: str) -> Dict:
    """Return the run's used_parameters.yml."""
    with open(os.path.join(run_directory(pattern), "used_parameters.yml")) as handle:
        return yaml.safe_load(handle)


def read_step_table(pattern: str) -> Dict:
    """Read timesteps.txt: step times, scale factors and update counts.

    ``Time-bins`` holds two columns in every row, which shifts every column
    after it by one; the columns are found from the table's own header.
    """
    path = os.path.join(run_directory(pattern), "timesteps.txt")
    with open(path) as handle:
        lines = handle.read().splitlines()
    header = [
        line.split()[1:]
        for line in lines
        if line.startswith("#") and "Time-bins" in line and "Updates" in line
    ][-1]
    i_step = header.index("Step")
    i_time = header.index("Time")
    i_a = header.index("Scale-factor")
    i_dt = header.index("Time-step")
    i_gas = header.index("Updates") + 1
    i_star = header.index("s-Updates") + 1
    rows = []
    for line in lines:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.split()
        rows.append(
            (
                int(fields[i_step]),
                float(fields[i_time]),
                float(fields[i_a]),
                float(fields[i_dt]),
                int(fields[i_gas]),
                int(fields[i_star]),
            )
        )
    table = np.array(rows, dtype=float)
    if not np.array_equal(table[:, 0], np.arange(len(table))):
        raise RuntimeError(f"{path}: the step numbers are not 0, 1, 2, ...")
    return {
        "path": path,
        "time": table[:, 1],
        "a": table[:, 2],
        "dt": table[:, 3],
        "gas": table[:, 4].astype(int),
        "star": table[:, 5].astype(int),
    }


def hubble_internal(a: np.ndarray, snap: Dict) -> np.ndarray:
    """Return H(a) in internal units (0 without cosmology)."""
    if not snap["cosmological"]:
        return np.zeros_like(np.asarray(a, dtype=float))
    c = snap["cosmology"]
    e2 = (
        c["Omega_r"] * a**-4
        + c["Omega_m"] * a**-3
        + c["Omega_k"] * a**-2
        + c["Omega_lambda"]
    )
    return c["H0 [internal units]"] * np.sqrt(e2)


def step_grid(run: List[Dict], parameters: Dict, table: Dict) -> Dict:
    """Return the uniform step grid of the run.

    With cosmology the time line is uniform in ln a and the steps here are a
    power-of-two subdivision of ``ln(a_end/a_begin)``; without it, of
    ``time_end - time_begin``. The subdivision is read from the table's first
    step and checked against every row.
    """
    first = run[0]
    if first["cosmological"]:
        span = np.log(
            float(parameters["Cosmology"]["a_end"])
            / float(parameters["Cosmology"]["a_begin"])
        )
        printed = np.diff(np.log(table["a"]))
    else:
        integration = parameters["TimeIntegration"]
        span = float(integration["time_end"]) - float(integration["time_begin"])
        printed = np.diff(table["time"])
    k = int(round(np.log2(span / printed[0])))
    step = span / 2**k
    # The table prints the scale factor to 7 decimals and the time to 7
    # significant digits: compare at that precision.
    if np.max(np.abs(printed - step)) > 0.05 * step:
        raise RuntimeError(
            f"{table['path']}: the steps are not one uniform subdivision of "
            f"the run span (first step {printed[0]:.6e}, grid {step:.6e})"
        )
    n_steps = len(table["time"]) - 1
    if first["cosmological"]:
        s = np.log(float(parameters["Cosmology"]["a_begin"])) + step * np.arange(
            n_steps + 1
        )
        # Friedmann time of every step boundary, from the start.
        mid = 0.5 * (s[1:] + s[:-1])
        nodes = mid[:, None] + 0.5 * step * GAUSS_NODES[None, :]
        dt = (
            0.5
            * step
            * np.sum(
                GAUSS_WEIGHTS[None, :] / hubble_internal(np.exp(nodes), first), axis=1
            )
        )
        t = np.concatenate(([0.0], np.cumsum(dt)))
    else:
        s = np.zeros(n_steps + 1)
        t = step * np.arange(n_steps + 1)
        dt = np.full(n_steps, step)
    return {"s": s, "t": t, "dt": dt, "step": step, "n_steps": n_steps}


def snapshot_steps(run: List[Dict], grid: Dict) -> Tuple[np.ndarray, float]:
    """Return the number of completed steps at each snapshot, and the clock gap.

    A snapshot on a step boundary may be written before or after that step's
    tasks, so the doses it holds are ambiguous: it gets -1 and is not compared.
    The exception is the last snapshot at the end of the run, written after the
    last step. The clock gap is the largest relative difference between SWIFT's
    elapsed time at a snapshot and the Friedmann time the models integrate.
    """
    first = run[0]
    completed = []
    clock = 0.0
    for snap in run:
        elapsed = (snap["time"] - first["time"]) / first["time_unit"]
        if first["cosmological"]:
            position = (np.log(snap["a"]) - grid["s"][0]) / grid["step"]
        else:
            position = elapsed / grid["step"]
        n = int(np.floor(position + 0.5))
        on_boundary = abs(position - n) < STEP_BOUNDARY_GUARD
        if snap is run[-1] and on_boundary and n == grid["n_steps"]:
            completed.append(n)
        elif snap is first:
            completed.append(0)
            continue
        elif on_boundary:
            print(
                f"  snapshot at a = {snap['a']:.9g} sits on a step boundary: "
                "skipped, the doses it holds are ambiguous"
            )
            completed.append(-1)
            continue
        else:
            n = int(np.floor(position))
            completed.append(n)
        if first["cosmological"] and n > 0:
            # Friedmann time at the snapshot, from the grid and one partial step.
            s0 = grid["s"][n]
            s1 = np.log(snap["a"])
            nodes = 0.5 * (s0 + s1) + 0.5 * (s1 - s0) * GAUSS_NODES
            partial = (
                0.5
                * (s1 - s0)
                * np.sum(GAUSS_WEIGHTS / hubble_internal(np.exp(nodes), first))
            )
            model = grid["t"][n] + partial
            clock = max(clock, abs(model / elapsed - 1.0))
    return np.array(completed), clock


def models(
    grid: Dict,
    run: List[Dict],
    eps: float,
    lam: Dict[str, float],
    lum: Dict[str, float],
    completed: np.ndarray,
    expansion: bool = True,
) -> Dict:
    """Return the discrete recursion, the continuous solution and E_noH.

    All three are per band at every snapshot, in internal units. The dose is
    ``L_b`` times each step narrowed to a float, as the star stores it.
    Without ``expansion`` both models drop their Hubble terms.
    """
    first = run[0]
    n_max = int(np.max(completed))
    dt = grid["dt"][:n_max]
    h_end = hubble_internal(np.exp(grid["s"][1 : n_max + 1]), first)
    if not expansion:
        h_end = np.zeros_like(h_end)
    dt32 = dt.astype(np.float32).astype(np.float64)
    # Dose folded at step n (1-based) is the star's dose of step n-1; the
    # fake step 0 injects nothing.
    folded = np.concatenate(([0.0], dt32[:-1]))

    rec = {b: np.zeros(n_max + 1) for b in BANDS}
    depth = {b: eps * lam[b] * h_end * dt32 for b in BANDS}
    for n in range(1, n_max + 1):
        i = n - 1
        lw_prev = rec["LW"][n - 1]
        for b in BANDS:
            a = depth[b][i]
            decay = np.exp(-a)
            phi = 1.0 if a == 0.0 else -np.expm1(-a) / a
            rec[b][n] = decay * rec[b][n - 1] + phi * eps * lum[b] * folded[i]
        a_lw = depth["LW"][i]
        if a_lw > 0.0:
            phi_lw = -np.expm1(-a_lw) / a_lw
            absorbed = lw_prev * (-np.expm1(-a_lw)) + eps * lum["LW"] * folded[i] * (
                1.0 - phi_lw
            )
            rec["PE"][n] += (lam["LW"] - 1.0) / lam["LW"] * absorbed

    no_h = {b: eps * lum[b] * np.concatenate(([0.0], np.cumsum(folded))) for b in BANDS}

    # Continuous solution: the deficit with respect to E_noH, as a fraction.
    # Emission over steps 1..n-1 (Friedmann time), decay from emission to the
    # end of step n.
    cont_fraction = {b: np.zeros(len(completed)) for b in BANDS}
    for k, n in enumerate(completed):
        if n < 2:
            continue
        s_lo = grid["s"][: n - 1]
        s_hi = grid["s"][1:n]
        width = s_hi - s_lo
        if first["cosmological"]:
            nodes = 0.5 * (s_lo + s_hi)[:, None] + 0.5 * width[:, None] * GAUSS_NODES
            jac = (
                0.5
                * width[:, None]
                * GAUSS_WEIGHTS
                / hubble_internal(np.exp(nodes), first)
            )
            x = grid["s"][n] - nodes
            if not expansion:
                x = np.zeros_like(x)
        else:
            jac = np.full((n - 1, len(GAUSS_NODES)), 0.0)
            jac[:] = 0.5 * grid["step"] * GAUSS_WEIGHTS
            x = np.zeros_like(jac)
        emitted = np.sum(jac)

        def loss(lam_b: float) -> float:
            return float(np.sum(-np.expm1(-eps * lam_b * x) * jac))

        cont_fraction["LW"][k] = loss(lam["LW"]) / emitted
        cont_fraction["N"][k] = loss(lam["N"]) / emitted
        transfer = 0.0
        if lam["LW"] != lam["PE"]:
            gain = float(
                np.sum(
                    (np.expm1(-eps * lam["PE"] * x) - np.expm1(-eps * lam["LW"] * x))
                    * jac
                )
            )
            transfer = (
                lum["LW"]
                / lum["PE"]
                * (lam["LW"] - 1.0)
                / (lam["LW"] - lam["PE"])
                * gain
            )
        cont_fraction["PE"][k] = (loss(lam["PE"]) - transfer) / emitted

    # Leading order of recursion minus continuous, relative to E_noH:
    # eps lambda d(ln a) [1 + q <x>], from the one-step fold-in delay and from
    # H held at the step end, with q = -(1/2) dlnH/dlna and <x> the dose-weighted
    # age in ln a. PE's lambda is its net coefficient.
    s_emit = 0.5 * (grid["s"][: n_max - 1] + grid["s"][1:n_max])
    age = np.zeros(len(completed))
    for k, n in enumerate(completed):
        if n >= 2:
            weights = folded[1:n]
            age[k] = np.sum(weights * (grid["s"][n] - s_emit[: n - 1])) / np.sum(
                weights
            )
    q = 0.0
    if first["cosmological"] and expansion:
        a_mid = np.exp(grid["s"][n_max // 2])
        q = (
            -0.5
            * a_mid
            * (
                np.log(hubble_internal(a_mid * 1.001, first))
                - np.log(hubble_internal(a_mid / 1.001, first))
            )
            / (a_mid * np.log(1.001**2))
        )
    net = {
        "LW": lam["LW"],
        "N": lam["N"],
        "PE": lam["PE"] - (lam["LW"] - 1.0) * lum["LW"] / lum["PE"],
    }
    step = grid["step"] if (first["cosmological"] and expansion) else 0.0

    out = {"rec": {}, "cont": {}, "no_h": {}, "estimate": {}}
    for b in BANDS:
        out["rec"][b] = rec[b][completed]
        out["no_h"][b] = no_h[b][completed]
        out["cont"][b] = out["no_h"][b] * (1.0 - cont_fraction[b])
        out["estimate"][b] = eps * net[b] * step * (1.0 + q * age)
    return out


def star_neighbours(pattern: str, run: List[Dict], kernel_gamma: float) -> int:
    """Return the largest number of gas particles inside the star's kernel radius."""
    files = sorted(glob.glob(pattern))
    times = [s["time"] for s in run]
    count = 0
    for name in files:
        with h5py.File(name, "r") as handle:
            time = float(np.atleast_1d(handle["/Header"].attrs["Time"])[0])
            h_star = float(handle["/PartType4/SmoothingLengths"][0])
        matches = [s for s, t in zip(run, times) if t == time * s["time_unit"]]
        if not matches:
            continue
        snap = matches[0]
        box = snap["box_comoving"]
        offset = snap["pos_comoving"] - snap["star_comoving"]
        offset -= box * np.round(offset / box)
        radius = kernel_gamma * h_star
        count = max(count, int(np.sum(np.sum(offset**2, axis=1) < radius**2)))
    return count


def float_bar(n_ngb: int) -> Dict[str, float]:
    """Return the float terms of the bar, relative, dimensionless."""
    n_round = max(n_ngb - 1, 0)
    growth = 1.0 - n_round * U32
    terms = {
        "step to float": U32,
        "dose to float": U32,
        "drain division": U32,
        "weight closure": np.inf if growth <= 0.0 else n_round * U32 / growth,
        "reconstruction": 4.0 * U32,
    }
    return terms


def lambda_print_term(lam: Dict[str, float], signal: Dict[str, float]) -> float:
    """Return the effect of the log's five-digit lambdas on the Hubble signal."""
    worst = 0.0
    for b in BANDS:
        relative = max(
            0.5 * 10 ** (np.floor(np.log10(abs(lam[k]))) - 4) / abs(lam[k])
            for k in (["PE", "LW"] if b == "PE" else [b])
        )
        worst = max(worst, relative * abs(signal[b]))
    return worst


def check(opt: argparse.Namespace) -> bool:
    """Run every premise and gate; return True on PASS."""
    run = load_run(opt.snapshots)
    first = run[0]
    parameters = read_parameters(opt.snapshots)
    table = read_step_table(opt.snapshots)
    ok = True

    def fail(message: str) -> None:
        nonlocal ok
        print(f"  FAIL: {message}")
        ok = False

    # Premises.
    if not read_isrf_propagation(opt.snapshots):
        fail("GEARFeedback:ISRF_propagation is off")
    if read_c_hyp_scheme(opt.snapshots) != "fixed_fraction":
        fail("the closed system needs ISRF_c_hyp_scheme: fixed_fraction")
    for snap in run:
        if snap.get("n_stars", 0) != 1:
            fail("the models are written for exactly one star")
            return False
        if np.any(snap["Z"] != 0.0):
            fail("the models need zero metallicity (kappa = 0)")
            return False
    lum = {"PE": first["L_PE"], "LW": first["L_LW"]}
    lum["N"] = lum["LW"]  # radiation_set_star_band_luminosity copies LW
    for snap in run:
        if snap["L_PE"] != lum["PE"] or snap["L_LW"] != lum["LW"]:
            fail(
                "the star's band luminosity changed during the run "
                f"(PE {lum['PE']!r} -> {snap['L_PE']!r}); the models need it constant"
            )
            return False
    if not (lum["PE"] > 0.0 and lum["LW"] > 0.0):
        fail("the star must emit in both bands")
        return False
    partial_gas = int(np.sum(table["gas"][1:] != len(first["mass"])))
    partial_star = int(np.sum(table["star"][1:] != 1))
    print(
        f"one-bin premise: {len(table['gas']) - 1} steps, {partial_gas} with a "
        f"partial gas update, {partial_star} without the star"
    )
    if partial_gas or partial_star:
        fail(
            "PREMISE VIOLATED: not every particle and the star are updated on "
            "every step, so the totals do not obey the closed system"
        )
        return False
    lam_log = read_band_edge_weights(opt.snapshots)
    if lam_log is None:
        fail("no band-edge weights in the run's log")
        return False
    lam = {"PE": lam_log["PE"], "LW": lam_log["LW"], "N": lam_log["N_LW"]}
    c_hyp = np.concatenate([s["c_hyp"] for s in run[1:]])
    if not np.all(c_hyp == c_hyp[0]) or not c_hyp[0] > 0.0:
        fail("the recorded propagation speed is not one positive value")
        return False
    eps = float(c_hyp[0]) / C_LIGHT_CGS
    for snap in run:
        for b in BANDS:
            if not np.all(np.isfinite(snap[U_KEYS[b]])):
                fail(f"non-finite {BAND_NAMES[b]} in a snapshot")
                return False
    expansion = not opt.no_expansion_in_models
    if first["cosmological"] and not expansion:
        print("  models WITHOUT the Hubble terms (wrong-input test)")

    grid = step_grid(run, parameters, table)
    completed, clock = snapshot_steps(run, grid)
    m = models(grid, run, eps, lam, lum, completed, expansion)
    m_true = m if expansion else models(grid, run, eps, lam, lum, completed, True)

    n_ngb = star_neighbours(opt.snapshots, run, opt.kernel_gamma)
    terms = float_bar(n_ngb)
    measured = {
        b: np.array(
            [
                np.sum(s["mass"] * s[U_KEYS[b]]) / (s["mass_unit"] * s["energy_unit"])
                for s in run
            ]
        )
        for b in BANDS
    }
    for b in BANDS:
        if measured[b][0] != 0.0:
            fail(f"the {BAND_NAMES[b]} is not zero at the start (seeded field)")
    use = completed >= 2
    if not np.any(use):
        fail("no snapshot holds a folded dose")
        return False
    last = int(np.flatnonzero(use)[-1])
    signal = {
        b: float(
            (m_true["no_h"][b][last] - m_true["cont"][b][last])
            / m_true["no_h"][b][last]
        )
        for b in BANDS
    }
    terms["lambda print"] = lambda_print_term(lam, signal)
    terms["clock"] = clock
    bar = float(sum(terms.values()))

    print(
        f"setup: cosmological={first['cosmological']}, a {first['a']:.6g} to "
        f"{run[-1]['a']:.6g}, {grid['n_steps']} steps, eps = c_pin/c = {eps:.6e}, "
        f"lambda_E(PE) {lam['PE']}, lambda_E(LW) {lam['LW']}, lambda_N(LW) "
        f"{lam['N']}, L_LW/L_PE {lum['LW'] / lum['PE']:.6f}, star neighbours {n_ngb}"
    )
    print(
        "float bar terms: "
        + ", ".join(f"{k} {v:.2e}" for k, v in terms.items())
        + f"; total {bar:.3e}"
    )

    for b in BANDS:
        no_h = m["no_h"][b][use]
        r1 = (measured[b][use] - m["rec"][b][use]) / no_h
        r2 = (measured[b][use] - m["cont"][b][use]) / no_h
        disc = np.abs(m["rec"][b][use] - m["cont"][b][use]) / no_h
        if not (np.all(np.isfinite(r1)) and np.all(np.isfinite(r2))):
            fail(f"{BAND_NAMES[b]}: non-finite residual")
            continue
        k1 = int(np.argmax(np.abs(r1)))
        k2 = int(np.argmax(np.abs(r2) / (bar + disc)))
        disc_last = float((m["rec"][b][last] - m["cont"][b][last]) / m["no_h"][b][last])
        print(
            f"{BAND_NAMES[b]}: Hubble signal at the end {signal[b]:+.4e}; "
            f"recursion minus continuous at the end {disc_last:+.3e}, leading-order "
            f"estimate {m['estimate'][b][last]:+.3e}"
        )
        if opt.table:
            print(
                "    a          step   measured-recursion  measured-continuous  signal"
            )
            for k in np.flatnonzero(use):
                e_noh = m["no_h"][b][k]
                print(
                    f"    {run[k]['a']:.7f} {completed[k]:6d} "
                    f"{(measured[b][k] - m['rec'][b][k]) / e_noh:+.3e}          "
                    f"{(measured[b][k] - m['cont'][b][k]) / e_noh:+.3e}           "
                    f"{(m_true['no_h'][b][k] - m_true['cont'][b][k]) / m_true['no_h'][b][k]:+.3e}"
                )
        line1 = f"  (S1) worst |measured - recursion|/E_noH {abs(r1[k1]):.3e}, bar {bar:.3e}"
        line2 = (
            f"  (S2) worst |measured - continuous|/E_noH {abs(r2[k2]):.3e}, "
            f"bar {bar + disc[k2]:.3e} (float {bar:.3e} + discretisation {disc[k2]:.3e})"
        )
        print(line1 + ("  PASS" if abs(r1[k1]) <= bar else "  FAIL"))
        print(line2 + ("  PASS" if abs(r2[k2]) <= bar + disc[k2] else "  FAIL"))
        if abs(r1[k1]) > bar:
            ok = False
        if abs(r2[k2]) > bar + disc[k2]:
            ok = False
        if first["cosmological"]:
            ratio = abs(signal[b]) / (bar + disc[-1])
            word = "required" if b != "PE" else "reported, not required"
            print(f"  signal over S2 bar {ratio:.3g} ({word})")
            if b != "PE" and not ratio > 1.0:
                fail(
                    f"{BAND_NAMES[b]}: the Hubble signal does not exceed the bar; "
                    "this leg cannot tell an expanding box from a static one"
                )
    return ok


def main() -> int:
    """Run the check and print the verdict."""
    opt = parse_options()
    ok = check(opt)
    print(f"RESULT: {'PASS' if ok else 'FAIL'}")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
