################################################################################
# This file is part of SWIFT.
# Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
Cross-run comparisons for ISRFShearAsymmetry. mode=gate: the self-calibrating
gate on A_centroid/A_spread/A_energy against a zero-shear control (2.8).
mode=pair: the M-S5 paired +v/-v grid comparison (contrast variant, report
only, 2.7).
"""

import argparse
import glob
import json
import os
import sys

import h5py
import numpy as np
from scipy.spatial import cKDTree

FLOORS = dict(A_centroid=0.01, A_spread=0.01, A_energy=0.01)
FACTOR = 3.0


def parse_options():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--mode", choices=["gate", "pair"], required=True)
    parser.add_argument("--runs", nargs="+", help="Gate mode: sheared run directories.")
    parser.add_argument("--zero", help="Gate mode: the matching zero-shear control run directory.")
    parser.add_argument("--plus", help="Pair mode: the +v_shear run directory.")
    parser.add_argument("--minus", help="Pair mode: the -v_shear run directory.")
    parser.add_argument("--grid-n", type=int, default=32)
    return parser.parse_args()


def load_metrics(run_dir):
    with open(os.path.join(run_dir, "shear_metrics.json")) as f:
        return json.load(f)


def mode_gate(opt):
    zero = load_metrics(opt.zero)
    print(f"=== Gate: zero-shear control = {opt.zero} ===")
    if zero.get("void"):
        raise RuntimeError(f"{opt.zero} (the zero-shear control) is itself VOID; "
                            "it cannot be used as the gate's calibration baseline.")
    overall_ok = True
    for run in opt.runs:
        m = load_metrics(run)
        print(f"\n-- {run} --")
        if m.get("moment_weighting") != zero.get("moment_weighting"):
            raise RuntimeError(
                f"moment_weighting mismatch: {run}={m.get('moment_weighting')!r} vs "
                f"{opt.zero}={zero.get('moment_weighting')!r} -- one of these "
                "shear_metrics.json files predates the negative-weight-clipping fix; "
                "rerun isrf_shear_asymmetry_check.py on both before gating."
            )
        if m.get("void"):
            print("  VOID (KH contamination or negative-weight share) -- "
                  "gate not meaningful for this run.")
            continue
        for band in ("FUV", "LW"):
            b = m["bands"][band]
            b0 = zero["bands"][band]
            for key in ("A_centroid", "A_spread", "A_energy"):
                a = abs(b[key])
                a0 = abs(b0[key])
                limit = max(FACTOR * a0, FLOORS[key])
                ok = a <= limit
                overall_ok &= ok
                print(f"  {band} {key}: |A|={a:.4e}  limit=max(3*{a0:.4e}, {FLOORS[key]})="
                      f"{limit:.4e}  -> {'PASS' if ok else 'FAIL'}")
    print(f"\nOverall gate: {'PASS' if overall_ok else 'FAIL'}")
    if not overall_ok:
        sys.exit(1)


def load_snapshot_last(run_dir):
    files = sorted(glob.glob(os.path.join(run_dir, "snap", "snapshot_*.hdf5")))
    with h5py.File(files[-1], "r") as f:
        pos = f["/PartType0/Coordinates"][:, :]
        h = f["/PartType0/SmoothingLengths"][:].astype(np.float64)
        u_fuv = f["/PartType0/FUVSpecificEnergies"][:].astype(np.float64)
        u_lw = f["/PartType0/LWSpecificEnergies"][:].astype(np.float64)
        boxsize = np.asarray(f["/Header"].attrs["BoxSize"], dtype=float).flatten()[0]
    return pos, h, u_fuv, u_lw, boxsize


def wendland_c2(q):
    """Wendland C2 kernel weight (3D), q = r/H, H = kernel support radius."""
    w = np.zeros_like(q)
    sel = q < 1.0
    t = 1.0 - q[sel]
    w[sel] = (21.0 / (2.0 * np.pi)) * t**4 * (4.0 * q[sel] + 1.0)
    return w


def sph_interpolate(pos, h, field, boxsize, grid_n):
    """Interpolate `field` onto a grid_n^3 periodic grid using each
    particle's own kernel_gamma*h support and a Wendland-C2 weight."""
    kernel_gamma = 1.9365
    tree = cKDTree(pos, boxsize=boxsize)
    coords = (np.arange(grid_n) + 0.5) / grid_n * boxsize
    gx, gy, gz = np.meshgrid(coords, coords, coords, indexing="ij")
    grid_pts = np.stack([gx.ravel(), gy.ravel(), gz.ravel()], axis=1)
    h_max = np.max(h) * kernel_gamma
    result = np.zeros(grid_pts.shape[0])
    weight_sum = np.zeros(grid_pts.shape[0])
    # Query all particles within h_max of each grid point (uniform search
    # radius, then re-weight per-particle with its own H = kernel_gamma*h).
    pairs = tree.query_ball_point(grid_pts, r=h_max, workers=-1)
    for i, neighbours in enumerate(pairs):
        if not neighbours:
            continue
        neighbours = np.array(neighbours)
        dx = grid_pts[i] - pos[neighbours]
        dx -= boxsize * np.round(dx / boxsize)
        r = np.sqrt(np.sum(dx**2, axis=1))
        H = kernel_gamma * h[neighbours]
        q = r / H
        w = wendland_c2(q) / H**3
        wsum = w.sum()
        if wsum > 0:
            result[i] = np.sum(w * field[neighbours]) / wsum
            weight_sum[i] = wsum
    return result.reshape(grid_n, grid_n, grid_n)


def mode_pair(opt):
    pos_p, h_p, fuv_p, lw_p, box_p = load_snapshot_last(opt.plus)
    pos_m, h_m, fuv_m, lw_m, box_m = load_snapshot_last(opt.minus)
    boxsize = box_p

    for band, field_p, field_m in (("FUV", fuv_p, fuv_m), ("LW", lw_p, lw_m)):
        grid_p = sph_interpolate(pos_p, h_p, field_p, boxsize, opt.grid_n)
        grid_m = sph_interpolate(pos_m, h_m, field_m, boxsize, opt.grid_n)
        # Map the -v run through x -> L - x (flip the first grid axis).
        grid_m_mapped = grid_m[::-1, :, :]
        u_max = np.max(np.abs(grid_p))
        A_pair = np.max(np.abs(grid_p - grid_m_mapped)) / u_max if u_max > 0 else float("nan")
        print(f"{band}: A_pair = {A_pair:.4e}  (report only, no threshold)")


def main():
    opt = parse_options()
    if opt.mode == "gate":
        mode_gate(opt)
    else:
        mode_pair(opt)


if __name__ == "__main__":
    main()
