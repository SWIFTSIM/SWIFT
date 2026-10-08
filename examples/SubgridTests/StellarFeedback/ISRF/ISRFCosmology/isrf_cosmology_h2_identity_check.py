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
"""Compare the shielded H2 fixture run without cosmology and at high redshift.

The shielded photodissociation of ``h2_shielded`` (a star in dust-free gas,
Grackle's local H2 self-shielding) has no closed form, so the reference is the
same physical problem run without cosmology. The cosmological runs start from
the same physical density, temperature, star and box, and last 25 yr, so
expansion changes the physical state by H t ~ 3e-8. Their results must match
the reference.

Compared, per particle, at the same elapsed proper time (linear interpolation
between snapshots), over the second half of the run, once the LW field has
crossed the box:

- the H2 destruction exponent D = ln[x_H2(t)/x_H2(0)];
- the physical LW specific energy;
- the temperature, the shielding length and the predicted shielding factor
  (the check of ``../ISRFH2Photodissociation``, Eqs. 2 to 5b).

Bars. The cosmological run differs from the reference by round-off: every
comoving quantity is converted with powers of a, and the proper-time steps
come from the cosmological time tables (they differ from the reference steps
by a few 1e-6). The reference rerun with another thread count
(``--replicate``) changes the round-off and the task order in the same way, so
the bar is twice that change, plus the interpolation error (measured on the
reference by predicting its odd snapshots from its even ones, a quarter of
which applies at half the spacing) and a float32 floor for derived quantities.
A run with halved steps is not a valid bar: the propagation speed is
``C h / dt``, so halving the step changes the physical solution.

Both runs must give the star the same age at the start (``star_age`` in
``run.sh``): a star of age 0 does not inject on the first step.

Self-shielding length diagnostic. With ``H2_self_shielding`` 3, Grackle
builds the column from the Jeans length. If Grackle used a comoving length at
redshift z, the column would be (1+z) times larger. The check prints the
destruction exponent ratio that would follow, so a failure can be told apart
from that case.
"""

import argparse
import glob
import os
import sys
from typing import Dict, List

import h5py
import numpy as np

sys.path.insert(
    0,
    os.path.join(
        os.path.dirname(os.path.abspath(__file__)), "..", "ISRFH2Photodissociation"
    ),
)
from isrf_h2_photodissociation_check import (  # noqa: E402
    GAMMA,
    K_BOLTZMANN_CGS,
    KERNEL_GAMMA_WENDLAND_C2,
    MU_METAL,
    M_H_CGS,
    jeans_shielding_length,
    shielding_factor,
)

# Snapshot fields are float32 and the cosmological run stores comoving values
# converted back with a^n: a derived quantity such as the Jeans length combines
# about ten rounded values, 10 x 6e-8.
FLOAT32_FLOOR = 1e-6

SPECIES = ["HI", "HII", "HeI", "HeII", "HeIII", "e", "HM", "H2I", "H2II"]


def parse_options() -> argparse.Namespace:
    """Parse the command line."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--reference", required=True, help="Non-cosmological snapshots")
    parser.add_argument(
        "--replicate", required=True, help="Same run, another thread count"
    )
    parser.add_argument(
        "--cosmo", required=True, nargs="+", help="Cosmological snapshot globs"
    )
    parser.add_argument("--h2-self-shielding", type=int, default=3, choices=[2, 3])
    return parser.parse_args()


def physical(dataset: h5py.Dataset, a: float, unit: float) -> np.ndarray:
    """Return a dataset in physical CGS."""
    exponent = float(np.atleast_1d(dataset.attrs["a-scale exponent"])[0])
    return dataset[:].astype(np.float64) * a**exponent * unit


def load_run(pattern: str, mode: int) -> Dict:
    """Read every snapshot of a run into per-snapshot arrays sorted by ID."""
    files = sorted(glob.glob(pattern))
    if len(files) < 3:
        raise RuntimeError(f"Need at least three snapshots for {pattern!r}")
    rows = []
    for name in files:
        with h5py.File(name, "r") as handle:
            units = handle["/Units"].attrs
            length = float(np.atleast_1d(units["Unit length in cgs (U_L)"])[0])
            mass = float(np.atleast_1d(units["Unit mass in cgs (U_M)"])[0])
            time = float(np.atleast_1d(units["Unit time in cgs (U_t)"])[0])
            a = float(np.atleast_1d(handle["/Header"].attrs["Scale-factor"])[0])
            cosmo = handle["/Cosmology"].attrs
            is_cosmo = int(np.atleast_1d(cosmo.get("Cosmological run", [0]))[0]) == 1
            gas = handle["/PartType0"]
            order = np.argsort(gas["ParticleIDs"][:])
            row = {
                "a": a,
                "elapsed": float(np.atleast_1d(handle["/Header"].attrs["Time"])[0])
                * time,
                "density": physical(gas["Densities"], a, mass / length**3)[order],
                "u": physical(gas["InternalEnergies"], a, (length / time) ** 2)[order],
                "u_LW": physical(gas["LWSpecificEnergies"], a, (length / time) ** 2)[
                    order
                ],
                "h": physical(gas["SmoothingLengths"], a, length)[order],
            }
            for species in SPECIES:
                row[species] = gas[species][:].astype(np.float64)[order]
            metals = gas["MetalMassFractions"][:].astype(np.float64)
            row["Z"] = (metals[:, -1] if metals.ndim == 2 else metals)[order]
        inverse_mu = (
            row["HI"]
            + row["HII"]
            + row["e"]
            + 0.25 * (row["HeI"] + row["HeII"] + row["HeIII"])
            + row["HM"]
            + 0.5 * (row["H2I"] + row["H2II"])
            + row["Z"] / MU_METAL
        )
        row["mu"] = 1.0 / inverse_mu
        row["T"] = (GAMMA - 1.0) * row["u"] * row["mu"] * M_H_CGS / K_BOLTZMANN_CGS
        if mode == 3:
            row["l_shield"] = jeans_shielding_length(
                row["T"], row["density"], row["mu"]
            )
        else:
            row["l_shield"] = KERNEL_GAMMA_WENDLAND_C2 * row["h"]
        n_h2 = row["H2I"] * row["density"] / (2.0 * M_H_CGS)
        row["column"] = 2.0 * n_h2 * row["l_shield"]
        n_total = row["density"] * inverse_mu / M_H_CGS
        row["f_shield"] = shielding_factor(row["column"], row["T"], n_total)
        row["n_total"] = n_total
        # SWIFT also dumps at time_end, which can repeat the last output time.
        if rows and row["elapsed"] == rows[-1]["elapsed"]:
            continue
        rows.append(row)
    start = rows[0]["elapsed"]
    for row in rows:
        row["elapsed"] -= start
    run = {
        "a": np.array([r["a"] for r in rows]),
        "t": np.array([r["elapsed"] for r in rows]),
    }
    for key in ["u_LW", "T", "l_shield", "f_shield", "column", "n_total"]:
        run[key] = np.array([r[key] for r in rows])
    run["D"] = np.array([np.log(r["H2I"] / rows[0]["H2I"]) for r in rows])
    return run


def interpolate(run: Dict, key: str, t: float) -> np.ndarray:
    """Interpolate a per-particle quantity linearly in elapsed time."""
    times = run["t"]
    j = int(np.clip(np.searchsorted(times, t) - 1, 0, len(times) - 2))
    w = (t - times[j]) / (times[j + 1] - times[j])
    return (1.0 - w) * run[key][j] + w * run[key][j + 1]


def metrics(run: Dict, ref: Dict, times: np.ndarray) -> Dict[str, float]:
    """Return the worst median relative difference of each quantity over times."""
    worst = {key: 0.0 for key in ["D", "u_LW", "T", "l_shield", "f_shield"]}
    for t in times:
        base = {key: interpolate(ref, key, t) for key in worst}
        other = {key: interpolate(run, key, t) for key in worst}
        arrived = base["u_LW"] > 1e-3 * np.max(base["u_LW"])
        destroyed = np.abs(base["D"]) > 0
        for key in worst:
            mask = (
                destroyed
                if key == "D"
                else (arrived if key == "u_LW" else np.ones_like(arrived))
            )
            if not np.any(mask):
                continue
            rel = np.median(np.abs(other[key][mask] / base[key][mask] - 1.0))
            # np.maximum, not the builtin: `max(0.0, nan)` returns 0.0, so a
            # non-finite value was dropped and the finiteness guard below
            # could never fire. The bars come from this same function, so a
            # non-finite reference shrank a bar instead of failing.
            worst[key] = float(np.maximum(worst[key], rel))
    return worst


def main() -> int:
    """Run the comparison."""
    opt = parse_options()
    ref = load_run(opt.reference, opt.h2_self_shielding)
    replicate = load_run(opt.replicate, opt.h2_self_shielding)
    cosmos = {
        pattern: load_run(pattern, opt.h2_self_shielding) for pattern in opt.cosmo
    }

    end = min(
        [ref["t"][-1], replicate["t"][-1]] + [c["t"][-1] for c in cosmos.values()]
    )
    # Second half of the run: the LW field has crossed the box (about
    # 4e8 s of 8e8 s at level 5), so the front, where relative differences
    # of a quantity rising from zero are unbounded, is not in the window.
    times = ref["t"][(ref["t"] >= 0.5 * end) & (ref["t"] <= end)]
    if np.max(ref["u_LW"][-1]) <= 0.0 or not np.all(np.abs(ref["D"][-1]) > 0.0):
        print("RESULT: FAIL (the LW field never destroys H2 in the reference)")
        return 1
    step_change = metrics(replicate, ref, times)

    half = {"t": ref["t"][::2]}
    for key in ["D", "u_LW", "T", "l_shield", "f_shield"]:
        half[key] = ref[key][::2]
    odd_times = ref["t"][1:-1:2]
    odd_times = odd_times[(odd_times >= 0.5 * end) & (odd_times <= end)]
    interpolation = metrics(half, ref, odd_times)

    print(
        f"h2_shielded identity, H2_self_shielding={opt.h2_self_shielding}, "
        f"{len(times)} comparison times up to {end:.4e} s"
    )
    print(
        f"  reference final median D {np.median(ref['D'][-1]):.4e}, f_shield "
        f"{np.median(ref['f_shield'][-1]):.4e}, l_shield "
        f"{np.median(ref['l_shield'][-1]) / 3.0856775814913673e18:.4g} pc"
    )
    bars = {}
    for key in step_change:
        bars[key] = 2.0 * step_change[key] + 0.25 * interpolation[key] + FLOAT32_FLOOR
        print(
            f"  bar {key:<9s}: 2 x replicate change {step_change[key]:.2e} + "
            f"interpolation {0.25 * interpolation[key]:.2e} + float32 {FLOAT32_FLOOR:.0e} "
            f"= {bars[key]:.2e}"
        )

    ok = True
    for pattern, run in cosmos.items():
        z = 1.0 / run["a"][0] - 1.0
        worst = metrics(run, ref, times)
        print(
            f"  run {pattern} (z = {z:.3g}, a span {run['a'][-1] / run['a'][0] - 1:.2e})"
        )
        for key in worst:
            passed = np.isfinite(worst[key]) and worst[key] <= bars[key]
            ok &= bool(passed)
            print(
                f"    {'PASS' if passed else 'FAIL'}: {key:<9s} worst median rel. diff "
                f"{worst[key]:.3e}, bar {bars[key]:.3e}"
            )
        if opt.h2_self_shielding == 3:
            t_end = times[-1]
            column = interpolate(ref, "column", t_end)
            temperature = interpolate(ref, "T", t_end)
            n_total = interpolate(ref, "n_total", t_end)
            f_phys = shielding_factor(column, temperature, n_total)
            f_comov = shielding_factor(column * (1.0 + z), temperature, n_total)
            print(
                f"    comoving Jeans length case: D would be {np.median(f_comov / f_phys):.3e} "
                "times the reference; measured ratio "
                f"{np.median(interpolate(run, 'D', t_end) / interpolate(ref, 'D', t_end)):.4f}"
            )
    print("RESULT: PASS" if ok else "RESULT: FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
