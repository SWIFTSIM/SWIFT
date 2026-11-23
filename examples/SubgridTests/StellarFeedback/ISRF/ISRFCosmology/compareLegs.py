#!/usr/bin/env python3
"""Print the energy budget of the light of the star, per band, for one or more runs.

Each entry is the sum over the particles of the last snapshot, as a fraction of the energy the star
has emitted: light still in flight (held, inbox and reservoir), absorbed by dust, redshifted by the
expansion, lost with no downstream neighbour, escaped from the star without reaching gas. The residual closes the
budget: it is zero up to float32 rounding, so the script exits with status 1 if |residual| exceeds 1e-6 in any
band of any run, or if any row is not finite. PE is the photoelectric band (6 to 11.2 eV), LW the Lyman-Werner band (11.2 to 13.6 eV).
The budget needs the Packet fields of a build configured with --enable-isrf-debugging-checks.
"""

import glob
import sys

import h5py
import numpy as np

RESIDUAL_BAR = 1e-6

GAS = (
    "PacketHeldEnergies",
    "PacketInboxEnergies",
    "PacketReservoirEnergies",
    "PacketAbsorbedEnergies",
    "PacketEscapedEnergies",
    "PacketRedshiftedEnergies",
)


def budget(run):
    """Return the scale factor and the budget rows of the last snapshot of a run."""
    files = sorted(glob.glob(f"{run}/snap/snapshot_*.hdf5"))
    if len(files) < 2:
        sys.exit(f"{run}: need at least two snapshots.")
    with h5py.File(files[-1], "r") as f:
        try:
            gas = {k: f["PartType0"][k][:].astype(np.float64).sum(axis=0) for k in GAS}
            star = f["PartType4"]
            emitted = star["PacketEmittedEnergies"][:].astype(np.float64).sum(axis=0)
            star_escaped = (
                star["PacketEscapedEnergies"][:].astype(np.float64).sum(axis=0)
            )
        except KeyError as err:
            sys.exit(
                f"{run}: field {err} missing: configure with --enable-isrf-debugging-checks."
            )
        a = (
            float(np.ravel(f["Cosmology"].attrs["Scale-factor"])[0])
            if "Cosmology" in f
            else 1.0
        )
    if not (np.all(np.isfinite(emitted)) and np.all(emitted > 0.0)):
        sys.exit(f"{run}: the star emitted no finite light.")
    inflight = (
        gas["PacketHeldEnergies"]
        + gas["PacketInboxEnergies"]
        + gas["PacketReservoirEnergies"]
    )
    rows = {
        "emitted (internal units)": emitted,
        "in flight": inflight / emitted,
        "absorbed": gas["PacketAbsorbedEnergies"] / emitted,
        "redshifted": gas["PacketRedshiftedEnergies"] / emitted,
        "no downstream neighbour": gas["PacketEscapedEnergies"] / emitted,
        "escaped from the star": star_escaped / emitted,
    }
    rows["residual"] = 1.0 - sum(
        v for k, v in rows.items() if k not in ("emitted (internal units)", "residual")
    )
    return a, rows


def main():
    runs = sys.argv[1:]
    if not runs:
        sys.exit("Usage: compareLegs.py RUN_DIR [RUN_DIR ...]")
    results = [budget(r) for r in runs]
    for band, name in enumerate(("PE", "LW")):
        print(f"{name} band")
        print(
            f"  {'':26s}"
            + "".join(f"{r.rstrip('/').split('/')[-1]:>16s}" for r in runs)
        )
        print(f"  {'scale factor':26s}" + "".join(f"{a:16.6f}" for a, _ in results))
        for key in results[0][1]:
            print(
                f"  {key:26s}"
                + "".join(f"{rows[key][band]:16.4e}" for _, rows in results)
            )
    failed = False
    for r, (_, rows) in zip(runs, results):
        for key, v in rows.items():
            if not np.all(np.isfinite(v)):
                print(f"FAIL {r}: non-finite value in '{key}'")
                failed = True
        if not np.all(np.abs(rows["residual"]) <= RESIDUAL_BAR):
            print(f"FAIL {r}: |residual| above {RESIDUAL_BAR:.0e}")
            failed = True
    if failed:
        sys.exit(1)


if __name__ == "__main__":
    main()
