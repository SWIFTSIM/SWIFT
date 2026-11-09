#!/usr/bin/env python3
"""Check the depth of the shadow that a clump casts in the ISRF.

Usage: checkShadow.py [CLUMP_RUN [CONTROL_RUN]]   (defaults: ., ../run_clump0)

CLUMP_RUN and CONTROL_RUN are the run directories of clump=1 and clump=0. Both hold
geometry.json and snap/snapshot_*.hdf5. The result is read from the last snapshot.

Kept depth, per zone: (1 - ratio) / (1 - expected).
  ratio     sum of the PE-band field with the clump over the same sum without it,
            over the gas of the zone.
  expected  the same ratio for an exact shadow: the analytic transmission of the clump
            along the ray from each star, averaged with the no-clump field of each
            particle as weight, the weight that the measured ratio carries.
An exact shadow scores 1. The bar is 0.90: the run passes when the shadow keeps at
least 90 per cent of the analytic depth. The bar is one-sided, so a deeper shadow
passes. It is a fixed number, not a fit to any run.

The zones lie behind the clump, between R + skirt and the zone radius from its centre.
A gas particle is in the zone of a star when its ray from that star crosses the core
of the clump (transmission below 0.01) and, with two stars, when the segment from the
other star stays at least R + skirt from the clump centre: a segment that grazes the
smoothed edge of the clump is partly absorbed. The zone of both stars (overlap) needs
a ray from each star through the core.

The reference field is the exact free-streaming field of the stars, with the ambient
extinction and no image of any star. The no-clump run is compared with it: the ratio
is 1 once the light front has passed and the field is steady. A ratio outside
[0.95, 1.05] prints a warning. It is not a gate.

The check fails (exit code 2) when a value is not finite, when a zone holds fewer than
100 particles, when the field of a zone has not settled (the no-clump field is below
0.6 of free streaming: the light front has not arrived) or when a kept depth is below
the bar.
"""

import glob
import json
import math
import os
import sys

import h5py
import numpy as np

BAR = 0.90
MIN_ZONE = 100
SETTLED = 0.6
CONTROL_BAND = (0.95, 1.05)
CORE_TRANSMISSION = 0.01


def chord_and_impact(p, star, radius):
    """Return the chord of the segment star -> p inside the sphere at the origin, and the impact parameter of its line."""
    v = p - star
    length = np.sqrt((v * v).sum(1))
    t = -(star[None, :] * v).sum(1) / length**2
    closest = star[None, :] + t[:, None] * v
    b = np.sqrt((closest * closest).sum(1))
    half = np.sqrt(np.maximum(radius**2 - b**2, 0.0))
    chord = np.clip(
        np.minimum(t * length + half, length) - np.maximum(t * length - half, 0.0),
        0.0,
        None,
    )
    return np.where(b < radius, chord, 0.0), b


def segment_distance(p, star):
    """Return the distance from the origin to the segment star -> p."""
    v = p - star
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.clip(-(star[None, :] * v).sum(1) / (v * v).sum(1), 0.0, 1.0)
    return np.linalg.norm(star[None, :] + t[:, None] * v, axis=1)


def zone_masks(rel, geo):
    """Return the zone masks and the analytic ratio of every particle.

    rel is the position relative to the clump centre, in h0. The ratio weights the
    stars by their ambient free-streaming field and applies the transmission of each.
    """
    radius = geo["clump_radius_h"]
    skirt = geo["skirt_h"]
    tau_per_h = geo["tau_diameter"] / (2.0 * radius)
    stars = (np.array(geo["stars_kpc"]) - np.array(geo["centre_kpc"])) / geo["h0_kpc"]
    transmission, weight, deep = [], [], []
    with np.errstate(divide="ignore", invalid="ignore"):
        for s in stars:
            chord, b = chord_and_impact(rel, s, radius)
            transmission.append(np.exp(-tau_per_h * chord))
            deep.append(
                (b <= radius - geo["edge_h"]) & (transmission[-1] < CORE_TRANSMISSION)
            )
            r = np.maximum(np.linalg.norm(rel - s, axis=1), 1e-9)
            weight.append(np.exp(-geo["kappa_ambient_per_h0"] * r) / r**2)
    r_c = np.linalg.norm(rel, axis=1)
    behind = (r_c >= radius + skirt) & (r_c <= geo["zone_radius_max_h"])
    expected = sum(w * t for w, t in zip(weight, transmission)) / sum(weight)
    if len(stars) == 1:
        return {"clump_shadow": deep[0] & behind}, expected
    clear = [segment_distance(rel, s) >= radius + skirt for s in stars]
    zones = {
        "star1_only": deep[0] & clear[1] & behind,
        "star2_only": deep[1] & clear[0] & behind,
        "overlap": deep[0] & deep[1] & behind,
    }
    return zones, expected


def free_streaming(x, stars, luminosity, c, rho, kappa):
    """Return the free-streaming specific energy L exp(-kappa r) / (4 pi r^2 c rho) of all stars at x."""
    u = np.zeros(len(x))
    for s, lum in zip(stars, luminosity):
        r = np.maximum(np.linalg.norm(x - s, axis=1), 1e-30)
        u += lum * np.exp(-kappa * r) / (4.0 * np.pi * r**2 * c * rho)
    return u


def read_snapshot(path):
    """Return the gas sorted by particle id, the stars and the time of a snapshot."""
    with h5py.File(path, "r") as f:
        gas, stars = f["PartType0"], f["PartType4"]
        order = np.argsort(gas["ParticleIDs"][:], kind="stable")
        u = np.asarray(gas["PESpecificEnergies"][:], dtype=np.float64)
        return dict(
            ids=np.asarray(gas["ParticleIDs"][:])[order],
            x=np.asarray(gas["Coordinates"][:], dtype=np.float64)[order],
            u=u.reshape(len(u), -1)[:, 0][order],
            stars=np.asarray(stars["Coordinates"][:], dtype=np.float64),
            luminosity=np.asarray(stars["PELuminosities"][:], dtype=np.float64).reshape(
                -1
            ),
            time=float(np.ravel(f["Header"].attrs["Time"])[0]),
            c=float(
                np.ravel(f["PhysicalConstants/InternalUnits"].attrs["speed_light_c"])[0]
            ),
        )


def gas_in_both(a, b):
    """Return both snapshots reduced to the gas present in both: a particle that left through a face is absent from one run."""
    ids, ia, ib = np.intersect1d(
        a["ids"], b["ids"], assume_unique=True, return_indices=True
    )
    if len(ids) < 0.9 * min(len(a["ids"]), len(b["ids"])):
        return None, None
    a, b = dict(a), dict(b)
    for d, idx in ((a, ia), (b, ib)):
        for key in ("ids", "x", "u"):
            d[key] = d[key][idx]
    return a, b


def zone_row(clump, control, geo):
    """Return the measured and analytic ratios of every zone of one snapshot pair."""
    rel = (clump["x"] - np.array(geo["centre_kpc"])) / geo["h0_kpc"]
    zones, ratio_exact = zone_masks(rel, geo)
    fs = free_streaming(
        control["x"],
        control["stars"],
        control["luminosity"],
        control["c"],
        geo["rho0_internal"],
        geo["kappa_ambient_per_kpc"],
    )
    rows = {}
    for name, mask in zones.items():
        n = int(mask.sum())
        row = dict(
            n=n,
            finite=bool(
                np.isfinite(fs[mask]).all() and np.isfinite(ratio_exact[mask]).all()
            ),
        )
        row["finite"] &= bool(
            np.isfinite(clump["u"][mask]).all()
            and np.isfinite(control["u"][mask]).all()
        )
        if n and row["finite"]:
            ctrl, cl, ref = (
                control["u"][mask].sum(),
                clump["u"][mask].sum(),
                fs[mask].sum(),
            )
            row["control_over_fs"] = ctrl / ref if ref > 0 else float("nan")
            row["ratio"] = cl / ctrl if ctrl > 0 else float("nan")
            row["expected"] = (
                float((control["u"][mask] * ratio_exact[mask]).sum() / ctrl)
                if ctrl > 0
                else float("nan")
            )
            with np.errstate(invalid="ignore", divide="ignore"):
                row["kept"] = (1.0 - row["ratio"]) / (1.0 - row["expected"])
            row["settled"] = bool(ref > 0 and ctrl >= SETTLED * ref)
        rows[name] = row
    return rows


def same_setup(clump_run, control_run):
    """Return the geometry of the clump run, or None after printing why the two runs differ."""
    geos = []
    for d in (clump_run, control_run):
        try:
            geos.append(json.load(open(os.path.join(d, "geometry.json"))))
        except OSError:
            print(f"FAIL: no geometry.json in {d}")
            return None
    for g in geos:
        g.pop("n_clump")
    if geos[0] != geos[1]:
        keys = [k for k in geos[0] if geos[0][k] != geos[1].get(k)]
        print(f"FAIL: geometry.json differs between the runs: {', '.join(keys)}")
        return None
    return geos[0]


def check(clump_run, control_run, passed=None):
    """Print the table of the kept depth and return the exit code: 0 pass, 2 fail.

    If passed is a dict, it receives the pass or fail of every zone.
    """
    geo = same_setup(clump_run, control_run)
    if geo is None:
        return 2
    files = [
        sorted(glob.glob(os.path.join(d, "snap", "snapshot_*.hdf5")))
        for d in (clump_run, control_run)
    ]
    if not files[0] or len(files[0]) != len(files[1]):
        print(
            f"FAIL: {len(files[0])} and {len(files[1])} snapshots: need the same number, at least one"
        )
        return 2
    for k, (fc, fn) in enumerate(zip(*files)):
        clump, control = read_snapshot(fc), read_snapshot(fn)
        last = k == len(files[0]) - 1
        if abs(clump["time"] - control["time"]) > 1e-6 * max(
            abs(control["time"]), 1e-30
        ):
            print(f"FAIL: {os.path.basename(fc)}: the runs are at different times")
            return 2
        clump, control = gas_in_both(clump, control)
        if clump is None:
            print(
                f"FAIL: {os.path.basename(fc)}: less than 90 per cent of the gas is in both runs"
            )
            return 2
        if not (
            clump["stars"].shape == control["stars"].shape
            and np.allclose(clump["stars"], control["stars"], rtol=1e-6, atol=0)
            and np.allclose(
                clump["luminosity"], control["luminosity"], rtol=1e-5, atol=0
            )
        ):
            print(
                f"FAIL: {os.path.basename(fc)}: the stars differ between the runs (position or PE luminosity)"
            )
            return 2
        if last and not (
            np.isfinite(clump["x"]).all() and np.isfinite(control["luminosity"]).all()
        ):
            print("FAIL: a position or a luminosity is not finite")
            return 2
        rows = zone_row(clump, control, geo)
        print(f"t = {clump['time']:.4e}, gas in both runs {len(clump['ids'])}")
        for name, row in rows.items():
            if "kept" in row:
                print(
                    f"  {name}: n = {row['n']}, ratio {row['ratio']:.4f}, analytic ratio {row['expected']:.4f}, "
                    f"kept {row['kept']:.4f}, control / free streaming {row['control_over_fs']:.3f}"
                )
            else:
                print(
                    f"  {name}: n = {row['n']}, no field"
                    + ("" if row["finite"] else " (not finite)")
                )
    failed = 0
    for name, row in rows.items():
        if not row["finite"]:
            verdict = "FAIL, a field is not finite"
        elif row["n"] < MIN_ZONE:
            verdict = f"FAIL, {row['n']} particles, fewer than {MIN_ZONE}"
        elif not row["settled"]:
            verdict = "FAIL, the light front has not arrived"
        elif not math.isfinite(row["kept"]):
            verdict = "FAIL, kept depth is not finite"
        else:
            verdict = f"{'pass' if row['kept'] >= BAR else 'FAIL'}, kept depth {row['kept']:.4f} against the bar {BAR}"
            if not CONTROL_BAND[0] <= row["control_over_fs"] <= CONTROL_BAND[1]:
                print(
                    f"WARNING {name}: control / free streaming {row['control_over_fs']:.3f} is outside "
                    f"[{CONTROL_BAND[0]}, {CONTROL_BAND[1]}]"
                )
        failed += verdict.startswith("FAIL")
        if passed is not None:
            passed[name] = verdict.startswith("pass")
        print(f"RESULT {name}: {verdict}")
    print("RESULT " + ("FAIL" if failed else "PASS"))
    return 2 if failed else 0


if __name__ == "__main__":
    args = sys.argv[1:]
    if len(args) > 2 or args[:1] in (["-h"], ["--help"]):
        print(__doc__)
        sys.exit(1)
    sys.exit(check(*(args + [".", "../run_clump0"][len(args) :])))
