#!/usr/bin/env python3
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
"""Check a run of the photoelectric equilibrium box against a one-zone integration.

The one-zone program (onezone.c, built against the same cooling library as swift) is
run at the median density, field and metallicity measured in the last snapshot. Its
equilibrium specific energy is the reference. The gates, with bars derived from the
mechanism and no fitted number:

  steady  the median specific energy changes over the last 20 per cent of the run by
          less than float32 rounding plus the work of the flow, per relaxation time
  equil   median specific energy against the one-zone value at the same time, within
          the decay left (times the starting offset), the work of the flow in one
          relaxation time and rounding, plus the float32 stall of the integrator
          (below)
  settled the one-zone value at the run's own end time is within the same bar of its
          value ten times later, so the run is long enough to compare an equilibrium;
          an unsettled run is excluded from the density ladder
  field   the seeded field is nonzero and identical in every snapshot (double precision)
  control (optional) a run without field ends at least ten times colder (specific
          energy, proportional to temperature)

Float32 stall. The integrator adds the cooling update to the specific energy in float32.
Each step moves u a fraction of about dt/tau of the way to its equilibrium. An update
smaller than half a unit in the last place (ulp) of u is rounded away, so u stops
changing once |u/u_eq - 1| * dt/tau < 0.5 ulp(u)/u. The gap to the one-zone value can
therefore be as large as 0.5 (ulp(u)/u) tau/dt, with dt the step of the run and tau the
one-zone relaxation time. The bar carries this term. It follows from the update rule and
is not fitted: it shrinks when the step is not small against tau, and halving the step
doubles it.

The output also gives the per-particle spread of the specific energy at the end. A ratio
of the 95th to the 5th percentile above 1.5 means the box split into two populations.
"""

import argparse
import glob
import os
import subprocess
import sys

import h5py
import numpy as np

SPEED_OF_LIGHT_CGS = 2.99792458e10
HABING_FLUX_CGS = 1.6e-3
PROTON_MASS_CGS = 1.67262192369e-24
HYDROGEN_MASS_FRACTION = 0.76
MYR_CGS = 3.15576e13
FLOAT32_EPS = 2.0**-24
GAMMA = 5.0 / 3.0


def parse_options():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("run_dir", help="Directory holding snap/ of the run")
    parser.add_argument("--onezone", required=True, help="Built onezone program")
    parser.add_argument("--table", required=True, help="Cloudy table file")
    parser.add_argument(
        "--network", type=int, default=1, help="Species network of the run"
    )
    parser.add_argument(
        "--pe-model",
        type=int,
        default=2,
        help="Library heating model number (2 is constant efficiency)",
    )
    parser.add_argument("--initial-temperature", type=float, default=1000.0)
    parser.add_argument(
        "--control", default=None, help="Directory of the run without field"
    )
    parser.add_argument(
        "--constant-g0",
        type=float,
        default=None,
        help="G0 of a run with a constant library background instead of a field",
    )
    parser.add_argument(
        "--clean-library-path",
        action="store_true",
        help="Remove LD_LIBRARY_PATH before running onezone (a laptop with an old Grackle on it)",
    )
    parser.add_argument(
        "--self-test", action="store_true", help="Feed the gates inputs that must fail"
    )
    return parser.parse_args()


def read_run(run_dir):
    """Return times [Myr], per-snapshot fields and the unit factors of a run."""
    files = sorted(glob.glob(os.path.join(run_dir, "snap", "*.hdf5")))
    if len(files) < 6:
        raise RuntimeError(f"Fewer than 6 snapshots in {run_dir}/snap")
    snaps = []
    for name in files:
        with h5py.File(name, "r") as f:
            units = f["Units"].attrs
            u_mass = float(np.ravel(units["Unit mass in cgs (U_M)"])[0])
            u_len = float(np.ravel(units["Unit length in cgs (U_L)"])[0])
            u_time = float(np.ravel(units["Unit time in cgs (U_t)"])[0])
            g = f["PartType0"]
            snaps.append(
                dict(
                    t_myr=float(np.ravel(f["Header"].attrs["Time"])[0])
                    * u_time
                    / MYR_CGS,
                    u=g["InternalEnergies"][:].astype(np.float64),
                    rho=g["Densities"][:].astype(np.float64),
                    pe=(
                        g["PESpecificEnergies"][:].astype(np.float64)
                        if "PESpecificEnergies" in g
                        else np.zeros(g["Densities"].shape)
                    ),
                    lw=(
                        g["LWSpecificEnergies"][:].astype(np.float64)
                        if "LWSpecificEnergies" in g
                        else np.zeros(g["Densities"].shape)
                    ),
                    div=g["VelocityDivergences"][:].astype(np.float64),
                    z=g["MetalMassFractions"][:, -1].astype(np.float64),
                )
            )
    conv = dict(
        density=u_mass / u_len**3,
        energy_per_mass=(u_len / u_time) ** 2,
        time=u_time,
    )
    return snaps, conv


def read_step_myr(run_dir, conv, t_from_myr):
    """Return the median time step [Myr] of the steps after t_from_myr, from timesteps.txt."""
    path = os.path.join(run_dir, "timesteps.txt")
    if not os.path.isfile(path):
        raise RuntimeError(f"{path} is missing: the float32 stall bar needs the step")
    data = np.loadtxt(path, usecols=(1, 4), comments="#", ndmin=2)
    finite_or_die("timesteps", [data])
    t_myr = data[:, 0] * conv["time"] / MYR_CGS
    dt_myr = data[:, 1] * conv["time"] / MYR_CGS
    keep = (t_myr >= t_from_myr) & (dt_myr > 0.0)
    if not np.any(keep):
        raise RuntimeError(f"No step with dt > 0 in {path} after {t_from_myr} Myr")
    return float(np.median(dt_myr[keep]))


def run_onezone(opt, n_h, g0, z, t_end_myr):
    cmd = [
        opt.onezone,
        str(opt.network),
        str(opt.pe_model),
        repr(n_h),
        repr(g0),
        repr(z),
        repr(opt.initial_temperature),
        "1e-20",
        repr(t_end_myr),
        opt.table,
        repr(0.01 / n_h),
    ]
    env = dict(os.environ)
    if opt.clean_library_path:
        env.pop("LD_LIBRARY_PATH", None)
    out = subprocess.run(
        cmd, capture_output=True, text=True, env=env, check=True
    ).stdout
    last = [line for line in out.splitlines() if line.startswith("FINAL")][-1]
    return {k: float(v) for k, v in (item.split("=") for item in last.split()[1:])}


def finite_or_die(name, arrays):
    for a in arrays:
        if not np.all(np.isfinite(a)):
            print(f"FAIL non-finite value in {name}")
            sys.exit(1)


def evaluate(snaps, conv, ref, label, dt_myr, ref_long=None, constant_field=False):
    """Return a list of (gate, value, bar, passed) for one run against the reference."""
    last = snaps[-1]
    tau_myr = ref["tau_Myr"]
    t_end = last["t_myr"]
    n_win = max(2, int(round(0.2 * len(snaps))))
    win = snaps[-n_win:]
    dt_win = win[-1]["t_myr"] - win[0]["t_myr"]

    u0 = np.median(win[0]["u"])
    u1 = np.median(win[-1]["u"])
    steady = abs(np.log(u1 / u0)) * tau_myr / dt_win
    div_myr = np.median(np.abs(last["div"])) * MYR_CGS / conv["time"]
    b_pdv = (GAMMA - 1.0) * div_myr * tau_myr
    b_round = 4.0 * FLOAT32_EPS
    b_steady = b_round * tau_myr / dt_win + b_pdv

    u_cgs = np.median(last["u"]) * conv["energy_per_mass"]
    equil = abs(u_cgs / ref["u_cgs"] - 1.0)
    offset = abs(
        np.median(snaps[0]["u"]) * conv["energy_per_mass"] / ref["u_cgs"] - 1.0
    )
    b_conv = offset * np.exp(-t_end / tau_myr)
    b_equil_core = np.sqrt(b_conv**2 + b_pdv**2 + b_round**2)
    u_med = np.float32(np.median(last["u"]))
    b_stall = 0.5 * float(np.spacing(u_med)) / float(u_med) * tau_myr / dt_myr
    b_equil = b_equil_core + b_stall

    field1 = snaps[0]["pe"] + snaps[0]["lw"]
    field_change = max(
        np.max(np.abs(s["pe"] + s["lw"] - field1) / np.maximum(np.abs(field1), 1e-300))
        for s in snaps
    )
    if not np.all(field1 > 0.0):
        field_change = np.inf
    b_field = 4.0 * 2.0**-52

    gates = [
        (f"{label} steady", steady, b_steady, steady < b_steady),
        (f"{label} equilibrium", equil, b_equil, equil < b_equil),
    ]
    if ref_long is not None:
        settled = abs(ref["u_cgs"] / ref_long["u_cgs"] - 1.0)
        gates.append(
            (f"{label} settled", settled, b_equil_core, settled < b_equil_core)
        )
    if not constant_field:
        gates.append(
            (f"{label} field constant", field_change, b_field, field_change < b_field)
        )
    return gates, dict(
        b_pdv=b_pdv,
        b_conv=b_conv,
        b_equil=b_equil,
        b_equil_core=b_equil_core,
        b_stall=b_stall,
        u_cgs=u_cgs,
        div_myr=div_myr,
        equil=equil,
    )


def main():
    opt = parse_options()
    snaps, conv = read_run(opt.run_dir)
    finite_or_die(
        "snapshots", [s[k] for s in snaps for k in ("u", "rho", "pe", "lw", "div", "z")]
    )
    last = snaps[-1]
    dt_myr = read_step_myr(opt.run_dir, conv, 0.8 * last["t_myr"])
    finite_or_die("time step", [np.array([dt_myr])])

    rho_cgs = last["rho"] * conv["density"]
    n_h = float(np.median(rho_cgs) * HYDROGEN_MASS_FRACTION / PROTON_MASS_CGS)
    g0_all = (
        SPEED_OF_LIGHT_CGS
        * rho_cgs
        * (last["pe"] + last["lw"])
        * conv["energy_per_mass"]
        / HABING_FLUX_CGS
    )
    g0 = float(np.median(g0_all)) if opt.constant_g0 is None else opt.constant_g0
    z = float(np.median(last["z"]))
    t_end_myr = last["t_myr"]
    ref = run_onezone(opt, n_h, g0, z, t_end_myr)
    ref_long = run_onezone(opt, n_h, g0, z, max(10.0 * t_end_myr, 50.0 / n_h))
    finite_or_die(
        "one-zone", [np.array(list(ref.values())), np.array(list(ref_long.values()))]
    )

    print(f"run {opt.run_dir}")
    print(
        f"measured median n_H = {n_h:.6g} cm^-3, median G0 = {g0:.6g} Habing "
        f"(spread {'n/a' if g0 == 0 else format(np.std(g0_all) / g0, '.3g')}), Z = {z:.6g}"
    )
    print(
        f"one-zone: T_eq = {ref['T']:.6g} K, u_eq = {ref['u_cgs']:.6g} erg/g, tau = {ref['tau_Myr']:.4g} Myr"
    )
    print(
        f"run length {last['t_myr']:.4g} Myr = {last['t_myr'] / ref['tau_Myr']:.1f} tau"
    )

    constant = opt.constant_g0 is not None
    results, info = evaluate(snaps, conv, ref, "run", dt_myr, ref_long, constant)
    e_cgs = last["u"] * conv["energy_per_mass"]
    p5, p50, p95 = np.percentile(e_cgs, [5, 50, 95])
    print(
        f"per-particle u at the end: p5 {p5:.4g}, p50 {p50:.4g}, p95 {p95:.4g} erg/g, "
        f"p95/p5 = {p95 / p5:.4g}"
        + (" SPLIT into two populations" if p95 / p5 > 1.5 else "")
    )
    print(
        f"convergence of the reference: u at the run end / u at ten times later - 1 = "
        f"{ref['u_cgs'] / ref_long['u_cgs'] - 1.0:.3g}"
    )
    print(
        f"median u at the end = {info['u_cgs']:.6g} erg/g, median |div v| = {info['div_myr']:.3g} / Myr"
    )
    print(
        f"step dt = {dt_myr:.4g} Myr (tau/dt = {ref['tau_Myr'] / dt_myr:.4g}), float32 stall term "
        f"{info['b_stall']:.4g}, bar without it {info['b_equil_core']:.4g}"
    )

    if opt.control:
        csnaps, cconv = read_run(opt.control)
        cu = np.median(csnaps[-1]["u"]) * cconv["energy_per_mass"]
        ratio = cu / info["u_cgs"]
        results.append(("control colder", ratio, 0.1, ratio < 0.1))

    if opt.self_test:
        for factor, name in (
            (1.5, "u x 1.5"),
            (1.01, "u x 1.01"),
            (1.0 + 3.0 * info["b_equil"], "u x (1 + 3 bar)"),
        ):
            bad = [dict(s, u=s["u"] * factor) for s in snaps]
            r, _ = evaluate(bad, conv, ref, name, dt_myr, None, constant)
            results.append(
                (
                    f"self-test must fail: equilibrium, {name}",
                    r[1][1],
                    r[1][2],
                    not r[1][3],
                )
            )
        if not constant:
            b_field = 4.0 * 2.0**-52
            for factor, name in (
                (1.0 + 1e-3, "field x 1.001"),
                (1.0 + 3.0 * b_field, "field x (1 + 3 bar)"),
            ):
                bad = [dict(s) for s in snaps]
                bad[-1] = dict(bad[-1], pe=bad[-1]["pe"] * factor)
                r, _ = evaluate(bad, conv, ref, name, dt_myr)
                results.append(
                    (f"self-test must fail: {name}", r[2][1], r[2][2], not r[2][3])
                )

    failed = 0
    for name, value, bar, passed in results:
        print(
            f"{'PASS' if passed else 'FAIL'}  {name}: value {value:.4g}, bar {bar:.4g}"
        )
        failed += not passed
    sys.exit(1 if failed else 0)


if __name__ == "__main__":
    main()
