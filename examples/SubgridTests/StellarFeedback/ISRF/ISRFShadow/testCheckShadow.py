#!/usr/bin/env python3
"""Run checkShadow.py on synthetic fields built from the gas of an initial-conditions file.

Usage: testCheckShadow.py RUN_DIR

RUN_DIR holds ICs_isrf_shadow.hdf5 and geometry.json (any run directory after makeIC.py).
The synthetic field with the clump is the exact free-streaming field times the analytic
transmission of the clump, with a planted error in every case but the first. Each case
states the verdict that the check must give. The exit code is 0 when all of them do.
"""

import contextlib
import io
import json
import os
import shutil
import sys
import tempfile

import h5py
import numpy as np

import checkShadow as cs

C_LIGHT = 299792.458  # km/s, internal velocity units
LUMINOSITY = 1.0


def write_snapshot(root, gas_x, ids, u, stars):
    """Write a snapshot with the fields that checkShadow.py reads."""
    os.makedirs(os.path.join(root, "snap"))
    with h5py.File(os.path.join(root, "snap", "snapshot_0000.hdf5"), "w") as f:
        f.create_group("Header").attrs["Time"] = [1.5e-3]
        f.create_group("PhysicalConstants/InternalUnits").attrs["speed_light_c"] = [
            C_LIGHT
        ]
        g = f.create_group("PartType0")
        g["ParticleIDs"], g["Coordinates"], g["PESpecificEnergies"] = ids, gas_x, u
        s = f.create_group("PartType4")
        s["Coordinates"], s["PELuminosities"] = stars, np.full(len(stars), LUMINOSITY)


def main(run_dir):
    geo = json.load(open(os.path.join(run_dir, "geometry.json")))
    with h5py.File(os.path.join(run_dir, "ICs_isrf_shadow.hdf5"), "r") as f:
        x = np.asarray(f["PartType0/Coordinates"][:], dtype=np.float64)
        stars = np.asarray(f["PartType4/Coordinates"][:], dtype=np.float64)
    ids = np.arange(len(x))
    h0, centre = geo["h0_kpc"], np.array(geo["centre_kpc"])
    rel = (x - centre) / h0
    rel_stars = (stars - centre) / h0
    zones, _ = cs.zone_masks(rel, geo)
    tau_per_h = geo["tau_diameter"] / (2.0 * geo["clump_radius_h"])
    radius, skirt = geo["clump_radius_h"], geo["skirt_h"]
    per_star, transmission, impact = [], [], []
    for s, sr in zip(stars, rel_stars):
        per_star.append(
            cs.free_streaming(
                x,
                [s],
                [LUMINOSITY],
                C_LIGHT,
                geo["rho0_internal"],
                geo["kappa_ambient_per_kpc"],
            )
        )
        chord, b = cs.chord_and_impact(rel, sr, radius)
        transmission.append(np.exp(-tau_per_h * chord))
        impact.append(b)
    control = sum(per_star)

    def field(shift=0.0, grazing=1.0):
        """Return the field with the clump: transmission raised by shift on deep rays, grazing rays scaled."""
        total = 0.0
        for u, t, b in zip(per_star, transmission, impact):
            t = np.where(t < cs.CORE_TRANSMISSION, np.minimum(t + shift, 1.0), t)
            t = np.where((b >= radius) & (b < radius + skirt), t * grazing, t)
            total = total + u * t
        return total

    first = next(iter(zones.values())).nonzero()[0]
    nan_field = field()
    nan_field[first[0]] = np.nan
    keep = np.ones(len(x), dtype=bool)
    keep[first[50:]] = False
    exact = field(0.0)
    cases = [
        ("exact shadow", dict(clump=exact, control=control), 0),
        ("transmission raised by 0.05", dict(clump=field(0.05), control=control), 0),
        (
            "transmission raised by 0.20 (too shallow)",
            dict(clump=field(0.2), control=control),
            2,
        ),
        ("no shadow", dict(clump=control, control=control), 2),
        ("a NaN in the field of a zone", dict(clump=nan_field, control=control), 2),
        (
            "both fields at 0.2 of free streaming (light front not arrived)",
            dict(clump=0.2 * exact, control=0.2 * control),
            2,
        ),
        (
            "first zone holds 50 particles",
            dict(clump=exact, control=control, keep=keep),
            2,
        ),
    ]
    if len(stars) == 2:
        cases.append(
            (
                "shadow 0.15 too shallow, grazing rays of the other star dimmed",
                dict(clump=field(0.15, 0.7), control=control),
                2,
            )
        )
    wrong = 0
    for name, case, wanted in cases:
        mask = case.get("keep", np.ones(len(x), dtype=bool))
        tmp = tempfile.mkdtemp(dir=os.environ.get("TMPDIR"))
        try:
            for label in ("clump", "control"):
                write_snapshot(
                    os.path.join(tmp, label),
                    x[mask],
                    ids[mask],
                    case[label][mask],
                    stars,
                )
            shutil.copy(
                os.path.join(run_dir, "geometry.json"), os.path.join(tmp, "clump")
            )
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                code = cs.check(
                    os.path.join(tmp, "clump"), os.path.join(tmp, "control")
                )
        finally:
            shutil.rmtree(tmp)
        print(
            f"{name}: exit code {code}, wanted {wanted}"
            + ("" if code == wanted else "  <-- WRONG")
        )
        for line in out.getvalue().splitlines():
            if line.startswith("RESULT") or line.startswith("  "):
                print("    " + line.strip())
        wrong += code != wanted
    print(
        "all cases give the wanted verdict"
        if not wrong
        else f"{wrong} case(s) give the wrong verdict"
    )
    return 1 if wrong else 0


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print(__doc__)
        sys.exit(1)
    sys.exit(main(sys.argv[1]))
