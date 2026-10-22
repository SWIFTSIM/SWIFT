#!/usr/bin/env python3
"""
Closed energy ledger of one run, gated on the float32 rounding bound.

For each band, at every snapshot where energy was injected since the first,

    R = (E - E_0 + Abs - Abs_0 - Inj + Inj_0) / (Inj - Inj_0),
    E = sum m (u + pending),

where Abs and Inj are the cumulative absorbed and injected specific energies,
pending is the energy finer neighbours owe a particle on a longer step, and
the subscript 0 marks the first snapshot.
The pair operator conserves energy exactly, so R is float32 rounding only:
each particle updates u and its two accumulators at most once per step, which
bounds |R| by n_steps * 2^-23, n_steps read from timesteps.txt.

The plain sums close only for one propagation speed in the whole box
(ISRF_c_hyp_scheme: fixed_fraction), and Abs and Inj are written only by a
build with --enable-debugging-checks. The check stops on either, and on any
missing or non-finite field, before it compares anything.

Exit status: 0 within the bound, 1 invalid input, 2 bound exceeded.
"""

import argparse
import glob
import sys

import h5py
import numpy as np

BANDS = ("PE", "LW")


def parse_options() -> argparse.Namespace:
    """Parse the command line."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("-s", "--snapshot", default="snap/snapshot_*.hdf5")
    parser.add_argument("--timesteps-log", default="timesteps.txt")
    return parser.parse_args()


def count_steps(path: str) -> int:
    """Return the number of steps listed in a SWIFT timesteps file."""
    with open(path) as f:
        return sum(1 for line in f if line.strip() and not line.startswith("#"))


def band_sums(path: str) -> dict:
    """
    Return the mass-weighted sums of one snapshot, per band.

    Parameters
    ----------
    path : str
        Snapshot file.

    Returns
    -------
    dict
        Band name to (E, Abs, Inj), each a float64 sum over the gas.
    """
    out = {}
    with h5py.File(path, "r") as h:
        gas = h["PartType0"]
        m = gas["Masses"][:].astype(np.float64)
        c = gas["HyperbolicPropagationSpeeds"][:]
        if not (np.all(np.isfinite(m)) and np.all(np.isfinite(c))):
            sys.exit(f"Non-finite mass or speed in {path}.")
        if c.size and np.any(c != c[0]):
            sys.exit(
                f"{path}: the propagation speed is not uniform; the plain sums "
                "close only with ISRF_c_hyp_scheme: fixed_fraction."
            )
        for band in BANDS:
            fields = []
            for name in (
                "SpecificEnergies",
                "PendingSpecificEnergies",
                "CumulativeAbsorbedSpecificEnergies",
                "CumulativeInjectedSpecificEnergies",
            ):
                key = band + name
                if key not in gas:
                    sys.exit(f"{path}: missing field {key}.")
                x = gas[key][:].astype(np.float64)
                if not np.all(np.isfinite(x)):
                    sys.exit(f"{path}: non-finite values in {key}.")
                fields.append(x)
            u, pending, absorbed, injected = fields
            out[band] = (
                np.sum(m * (u + pending)),
                np.sum(m * absorbed),
                np.sum(m * injected),
            )
    return out


def main() -> int:
    """Compute the ledger per snapshot and gate it."""
    args = parse_options()
    files = sorted(glob.glob(args.snapshot))
    if len(files) < 2:
        print("Need at least two snapshots.")
        return 1
    n_steps = count_steps(args.timesteps_log)
    if n_steps < 1:
        print(f"No step in {args.timesteps_log}.")
        return 1
    bar = n_steps * 2.0**-23

    first = band_sums(files[0])
    worst = {band: 0.0 for band in BANDS}
    last = None
    print("snapshot  " + "  ".join(f"R_{band:<9s}" for band in BANDS))
    for path in files:
        sums = band_sums(path)
        row = []
        for band in BANDS:
            E0, absorbed0, injected0 = first[band]
            E, absorbed, injected = sums[band]
            d_injected = injected - injected0
            if d_injected > 0.0:
                r = (E - E0 + absorbed - absorbed0 - d_injected) / d_injected
                if not np.isfinite(r):
                    print(f"{path}: non-finite residual in {band}.")
                    return 1
                worst[band] = max(worst[band], abs(r))
                row.append(f"{r:+.3e}")
            else:
                row.append("n/a       ")
        last = sums
        print(f"{path[-9:-5]:8s}  " + "  ".join(row))

    if any(last[band][2] - first[band][2] <= 0.0 for band in BANDS):
        print(
            "No energy injected after the first snapshot: build with "
            "--enable-debugging-checks."
        )
        return 1

    status = 0
    for band in BANDS:
        ok = worst[band] <= bar
        status = status if ok else 2
        print(
            f"{band}: max |R| {worst[band]:.3e}, bound {bar:.3e} "
            f"({n_steps} steps): {'pass' if ok else 'FAIL'}"
        )
    return status


if __name__ == "__main__":
    sys.exit(main())
