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
"""
Check that the ISRF dose-reservoir total does not depend on the injection
cadence.

The sweep is three runs of the same fixture at the same end time, whose
only difference is dt_max: cadence_base, cadence_half and cadence_quarter
put every gas particle and the star on one shared bin of, respectively,
dt, dt/2 and dt/4, so the number of injection events over the run changes
by a factor four while nothing else does.

With Z = 0 and a pinned c_hyp the stored energy E = (c / c_hyp) sum m u
trails the emitted energy L t by the dose still held between the star and
the gas. The gated quantity is that lag EXPRESSED IN THE RUN'S OWN STEP,
n = (L t - E) / (L dt). The scheme's hand-off pipeline is a fixed number
of steps deep, so n is the cadence-invariant quantity: the same band of
integers must come out of every run, whatever dt is.

That is what makes this check different from the bounded-window gate of
isrf_timestep_hierarchy_injection_check.py, which bounds n inside [-1, 3]
in each run separately. Three defect models are separated here:

  (a) dose legitimately in flight when the run stops: the lag scales with
      dt, so n is constant. This is what the runs should show.
  (b) a fixed fraction alpha of every lump dropped: the lag is alpha L t,
      cadence-independent, so n grows as alpha t / dt and the runs
      disagree by a factor four across the sweep. A per-run window cannot
      see this below its own 3-step bound; the cross-cadence comparison
      bounds alpha about a hundred times tighter.
  (c) a fixed absolute amount lost per injection event: the lag grows as
      1/dt and n as 1/dt^2, with the finest run worst.

Every gated value is tested for finiteness before it is compared, and the
script exits nonzero on any failure.
"""

import argparse
import glob
import os
import re
import sys

import h5py
import numpy as np

C_LIGHT_CGS = 2.99792458e10

# Relative rounding of the float32 specific energies the snapshots carry.
FLT_EPSILON = float(np.finfo(np.float32).eps)

# Pipeline depth of the dose hand-off, in steps of a uniform cadence. The
# emitted energy reaches the gas through three stages, each holding at
# most one cadence of dose: emission awaiting the star's next hand-off,
# the lump awaiting the receiving particle's next feedback_reset_part, and
# the drawn dose awaiting its release window. A lump is handed out for a
# star step that has not elapsed, which is a lead of one step, so the band
# is [-1, 3]. This is lag_window() of the companion check with every
# cadence set equal; see its docstring for the code references.
N_LAG_MIN = -1.0
N_LAG_MAX = 3.0

# Sweep members, coarsest first. The first is the reference.
RUNS = ("cadence_base", "cadence_half", "cadence_quarter")

# used_parameters.yml keys the sweep is allowed to differ in.
SWEPT_KEYS = ("dt_max",)


def parse_options() -> argparse.Namespace:
    """Parse the command line."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--dir", default=".", help="Directory holding the runs")
    return parser.parse_args()


def snapshots(run: str) -> list:
    """Return the sorted snapshot files of a run, one per output time.

    Parameters
    ----------
    run : str
        Run directory.

    Returns
    -------
    list
        Snapshot paths, in time order.
    """
    files, times = [], set()
    for path in sorted(glob.glob(os.path.join(run, "snap", "snapshot_*.hdf5"))):
        with h5py.File(path, "r") as f:
            t = float(np.asarray(f["/Header"].attrs["Time"]).flat[0])
        if not any(abs(t - s) <= 1e-9 * max(abs(t), 1e-300) for s in times):
            times.add(t)
            files.append(path)
    if len(files) < 3:
        raise RuntimeError(f"Fewer than three snapshots in {run}/snap")
    return files


def load(path: str) -> dict:
    """Read one snapshot.

    Parameters
    ----------
    path : str
        Snapshot file.

    Returns
    -------
    dict
        Snapshot fields in internal units.
    """
    with h5py.File(path, "r") as f:
        units = f["/Units"].attrs
        gas = f["/PartType0"]
        star = f["/PartType4"]
        snap = dict(
            time=float(np.asarray(f["/Header"].attrs["Time"]).flat[0]),
            ul=float(np.asarray(units["Unit length in cgs (U_L)"]).flat[0]),
            ut=float(np.asarray(units["Unit time in cgs (U_t)"]).flat[0]),
            mass=gas["Masses"][:].astype(np.float64),
            Z=gas["MetalMassFractions"][:, -1].astype(np.float64),
            u={
                "PE": gas["PESpecificEnergies"][:].astype(np.float64),
                "LW": gas["LWSpecificEnergies"][:].astype(np.float64),
            },
            L={
                "PE": float(star["PELuminosities"][0]),
                "LW": float(star["LWLuminosities"][0]),
            },
        )
    snap["c"] = C_LIGHT_CGS * snap["ut"] / snap["ul"]
    return snap


def step_table(run: str) -> np.ndarray:
    """Return timesteps.txt rows: time, dt, gas updates, star updates.

    Parameters
    ----------
    run : str
        Run directory.

    Returns
    -------
    numpy.ndarray
        One row per step.
    """
    rows = []
    with open(os.path.join(run, "timesteps.txt")) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            v = line.split()
            rows.append([float(v[1]), float(v[4]), int(v[7]), int(v[9])])
    return np.array(rows)


def parameter(run: str, key: str) -> float:
    """Read a numerical parameter from used_parameters.yml.

    Parameters
    ----------
    run : str
        Run directory.
    key : str
        Parameter name.

    Returns
    -------
    float
        Its value.
    """
    with open(os.path.join(run, "used_parameters.yml")) as f:
        m = re.search(rf"{key}:\s*([0-9.eE+-]+)", f.read())
    return float(m.group(1))


def gate(name: str, value: float, bar: float) -> bool:
    """Print and return one upper-bound gate; non-finite values fail.

    Parameters
    ----------
    name : str
        Label to print.
    value : float
        Measured value.
    bar : float
        Upper bound, inclusive.

    Returns
    -------
    bool
        True if the value is finite and at or below a finite bar.
    """
    ok = bool(np.isfinite(value) and np.isfinite(bar) and value <= bar)
    print(f"  {name}: {value:.4e} <= {bar:.4e} -> {'PASS' if ok else 'FAIL'}")
    return ok


def gate_window(name: str, value: float, low: float, high: float) -> bool:
    """Print and return one two-sided gate; non-finite values fail.

    Parameters
    ----------
    name : str
        Label to print.
    value : float
        Measured value.
    low, high : float
        Window edges, inclusive.

    Returns
    -------
    bool
        True if the value is finite and inside a finite window.
    """
    finite = bool(np.isfinite(value) and np.isfinite(low) and np.isfinite(high))
    ok = bool(finite and low <= value <= high)
    print(
        f"  {name}: {value:.4e} in [{low:.4e}, {high:.4e}] -> "
        f"{'PASS' if ok else 'FAIL'}"
    )
    return ok


def uniform_cadence(run: str, n_gas: int) -> tuple:
    """Check that one run has a single shared time bin, and return it.

    The comparison across the sweep is only meaningful if each run's
    hand-off stages all run at the same cadence, so this requires every
    step to update every gas particle and the star, and every step to
    carry the same dt.

    Parameters
    ----------
    run : str
        Run directory.
    n_gas : int
        Number of gas particles.

    Returns
    -------
    tuple
        (ok, dt, number of steps).
    """
    table = step_table(run)[1:]
    dt = table[:, 1]
    spread = float(np.max(np.abs(dt - dt[0])) / max(abs(dt[0]), 1e-300))
    ok = bool(
        np.all(np.isfinite(dt))
        and dt[0] > 0.0
        and spread <= FLT_EPSILON
        and np.all(table[:, 2] == n_gas)
        and np.all(table[:, 3] > 0)
    )
    print(
        f"  {run}: {len(table)} steps, dt {dt[0]:.4e}, dt spread "
        f"{spread:.2e}, all-gas steps {int(np.sum(table[:, 2] == n_gas))}, "
        f"star steps {int(np.sum(table[:, 3] > 0))} -> "
        f"{'PASS' if ok else 'FAIL'}"
    )
    return ok, float(dt[0]), len(table)


def same_configuration(base: str, other: str) -> bool:
    """Check that two runs' used_parameters.yml differ only in the sweep.

    A cross-cadence comparison is invalidated by any other configuration
    difference, the dissipation parameters in particular: they enter the
    same sum m u the gate reads, so a stray override would be measured as
    an injection defect.

    Parameters
    ----------
    base, other : str
        Run directories.

    Returns
    -------
    bool
        True if every differing line belongs to a swept key.
    """
    lines = []
    for run in (base, other):
        with open(os.path.join(run, "used_parameters.yml")) as f:
            # The file opens with a run-time stamp comment, which differs
            # between any two runs and says nothing about their physics.
            lines.append(
                [ln.rstrip() for ln in f if not ln.lstrip().startswith("#")]
            )
    ok = len(lines[0]) == len(lines[1])
    offenders = []
    if ok:
        for a, b in zip(lines[0], lines[1]):
            if a == b:
                continue
            if not any(k in a for k in SWEPT_KEYS):
                offenders.append(a.strip())
                ok = False
    print(
        f"  {other} configuration against {base}: "
        f"{len(offenders)} unexpected differences"
        + (f" ({'; '.join(offenders[:3])})" if offenders else "")
        + f" -> {'PASS' if ok else 'FAIL'}"
    )
    return bool(ok)


def lag_in_steps(snaps: list, band: str, pin: float, dt: float) -> tuple:
    """Return the in-flight lag of one band, in units of the run's step.

    Parameters
    ----------
    snaps : list
        Snapshots of one run, first one dropped.
    band : str
        "PE" or "LW".
    pin : float
        Pinned c_hyp, internal units.
    dt : float
        The run's shared time step.

    Returns
    -------
    tuple
        (n, tol, L), with n the lag in steps at each snapshot, tol its
        float32 allowance in steps, and L the stellar luminosity series.
    """
    L = np.array([s["L"][band] for s in snaps])
    t = np.array([s["time"] for s in snaps])
    E = np.array([np.sum(s["mass"] * s["u"][band]) * s["c"] / pin for s in snaps])
    n = (t - E / L) / dt
    # E is a sum of N float32 specific energies, so its rounding error
    # grows as sqrt(N) eps E; carried to the lag through L that is
    # sqrt(N) eps t, and to n a further division by dt.
    n_gas = len(snaps[0]["mass"])
    tol = np.sqrt(n_gas) * FLT_EPSILON * t / dt
    return n, tol, L


def main() -> None:
    """Run every gate and exit nonzero on any failure."""
    opt = parse_options()
    path = lambda name: os.path.join(opt.dir, name)
    ok = True

    print("Uniform cadence, matched configuration and absorption-free setup")
    info = {}
    for name in RUNS:
        snap = load(snapshots(path(name))[-1])
        r = uniform_cadence(path(name), len(snap["mass"]))
        ok &= r[0]
        info[name] = dict(dt=r[1], n_steps=r[2])
        ok &= gate(f"{name} max metallicity", float(snap["Z"].max()), 0.0)
        if name != RUNS[0]:
            ok &= same_configuration(path(RUNS[0]), path(name))

    print("Lag in units of the run's own step")
    measured = {}
    for name in RUNS:
        snaps = [load(f) for f in snapshots(path(name))[1:]]
        pin = parameter(path(name), "ISRF_c_hyp_pin_for_debugging")
        dt = info[name]["dt"]
        t_end = snaps[-1]["time"]
        info[name]["t_end"] = t_end
        for b in ("PE", "LW"):
            n, tol, L = lag_in_steps(snaps, b, pin, dt)
            positive = bool(np.all(np.isfinite(L)) and np.all(L > 0.0))
            print(
                f"  {name} {b} luminosity finite and positive: "
                f"{'PASS' if positive else 'FAIL'}"
            )
            ok &= positive
            spread = float(np.max(np.abs(L - L[0])) / max(abs(L[0]), 1e-300))
            ok &= gate(f"{name} {b} luminosity spread", spread, FLT_EPSILON)
            # np.min/np.max, not a nan-ignoring variant: a non-finite lag
            # must reach the gate rather than be dropped from it.
            n_lo = float(np.min(n + tol))
            n_hi = float(np.max(n - tol))
            measured[(name, b)] = dict(lo=n_lo, hi=n_hi, tol=float(tol[-1]))
            ok &= gate_window(
                f"{name} {b} lag floor, steps", n_lo, N_LAG_MIN, N_LAG_MAX
            )
            ok &= gate_window(
                f"{name} {b} lag ceiling, steps", n_hi, N_LAG_MIN, N_LAG_MAX
            )
            print(
                f"  {name} {b}: dt {dt:.4e}, {len(n)} snapshots, lag "
                f"{n_lo:.4f} to {n_hi:.4f} steps, budget in flight at t_end "
                f"{n[-1] * dt / t_end:.3e} of the emitted total"
            )

    print("Cadence independence of the lag band")
    ref = RUNS[0]
    for name in RUNS[1:]:
        for b in ("PE", "LW"):
            a, c = measured[(ref, b)], measured[(name, b)]
            bar = a["tol"] + c["tol"]
            ok &= gate(
                f"{name} {b} lag ceiling against {ref}",
                abs(c["hi"] - a["hi"]),
                bar,
            )
            ok &= gate(
                f"{name} {b} lag floor against {ref}",
                abs(c["lo"] - a["lo"]),
                bar,
            )

    print("Reported, not gated: what the sweep bounds")
    ref_dt, ref_t = info[ref]["dt"], info[ref]["t_end"]
    fine = RUNS[-1]
    fine_dt = info[fine]["dt"]
    print(
        f"  cadence range {ref_dt / fine_dt:.1f}x, {info[ref]['n_steps']} to "
        f"{info[fine]['n_steps']} steps over the same end time {ref_t:.4e}"
    )
    for b in ("PE", "LW"):
        a, c = measured[(ref, b)], measured[(fine, b)]
        bar = a["tol"] + c["tol"]
        # Model (b), a fixed fraction alpha of every lump dropped: the lag
        # is alpha L t whatever dt is, so n differs between two cadences by
        # alpha t (1/dt_fine - 1/dt_ref). The gate above bounds that
        # difference by bar, which inverts to a bound on alpha.
        lever = ref_t * (1.0 / fine_dt - 1.0 / ref_dt)
        print(
            f"  {b}: dropped-fraction bound from the sweep "
            f"{bar / lever:.3e}, against {N_LAG_MAX * ref_dt / ref_t:.3e} "
            f"from the coarsest run's own window alone"
        )
        # Model (c), a fixed absolute loss per injection event: the finest
        # run's own window bounds it hardest, one event being L dt.
        print(
            f"  {b}: per-event loss bound {N_LAG_MAX * fine_dt / ref_t:.3e} "
            f"of the emitted total, {N_LAG_MAX:.0f} lumps of the finest "
            f"cadence"
        )

    print("OVERALL:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
