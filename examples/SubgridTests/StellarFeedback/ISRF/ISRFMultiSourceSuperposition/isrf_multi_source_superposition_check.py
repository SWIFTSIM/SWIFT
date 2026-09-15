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
Check the ISRF injection and propagation with several sources.

Gates (each per band, FUV and LW):

E   Exact injection superposition, propagation off: per gas particle,
    |u(AB) - u(A) - u(B)| <= 1e-5 (u(A) + u(B)) + 4 eps (max u(A) +
    max u(B)), eps the float32 machine epsilon. The first term is float32
    storage; the second is the float32 evaluation of the kernel polynomial,
    whose absolute error near the support radius is a few eps of the
    central value (clamped at zero when negative), so a particle at
    r/H > 0.98 can receive a noise-level deposit in one run and none in
    another.
G1  Amplitude identity, propagation on, last snapshot:
    e = (sum_i m_i u_i / lambda_i - P / c) / (P / c), with
    P = sum_stars sum_j w_j L ext_j. Single sources: |e| <= 0.10.
    Pair: |e_AB - (e_A + e_B)/2| <= |e_A - e_B|/2 + 3 sigma, sigma the
    largest standard deviation of e over the last four snapshots of A, B.
L1  Amplitude identity of the source lattice: |e_lat| <= |e_single| +
    3 sigma_single (the one-star run at the same metallicity).
L2  Isotropic closure branch in the lattice: median reduced flux
    f = |F| / (c_hyp u) over gas outside every star's kernel <= 0.18,
    where the M1 Eddington factor is 5% above the isotropic value 1/3.

Reported without a bar: u(AB) - [u(A) + u(B)] with propagation on (the
M1 closure merges crossing beams, so the propagated field is not
additive), the lattice's spatial scatter of u against the continuum
free-streaming and diffusive lattice sums, and the one-star run's f.
Every gated value must be finite, or the gate fails.
"""

import argparse
import glob
import os
import re
import sys

import h5py
import numpy as np
from scipy.spatial import cKDTree

# src/feedback/GEAR/radiation.h
SIGMA_D_CGS = {"FUV": 9e-22, "LW": 1.5e-21}
MU_H = 1.4
M_H_CGS = 1.6726219e-24
GRACKLE_SOLAR_Z = 0.01295
# Wendland C2 kernel support over smoothing length, 3D
GAMMA_3D = 1.936492
C_LIGHT_CGS = 2.99792458e10
F_ISOTROPIC_BAR = 0.18
# Below this |F|, F.F underflows in float32 and the closure turns isotropic.
F_UNDERFLOW = 1.08e-19


def parse_options() -> argparse.Namespace:
    """Parse the command line."""
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--dir", default=".", help="Directory holding the runs")
    parser.add_argument("--window", type=int, default=4, help="Steady snapshots")
    parser.add_argument(
        "--sample", type=int, default=400, help="Particles for lattice sums"
    )
    return parser.parse_args()


def load(path: str) -> dict:
    """Read one snapshot, gas sorted by particle ID.

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
        order = np.argsort(gas["ParticleIDs"][:])
        star = f["/PartType4"]
        snap = dict(
            time=float(np.asarray(f["/Header"].attrs["Time"]).flat[0]),
            boxsize=float(np.asarray(f["/Header"].attrs["BoxSize"]).flat[0]),
            ul=float(np.asarray(units["Unit length in cgs (U_L)"]).flat[0]),
            um=float(np.asarray(units["Unit mass in cgs (U_M)"]).flat[0]),
            ut=float(np.asarray(units["Unit time in cgs (U_t)"]).flat[0]),
            pos=gas["Coordinates"][:][order].astype(np.float64),
            rho=gas["Densities"][:][order].astype(np.float64),
            h=gas["SmoothingLengths"][:][order].astype(np.float64),
            mass=gas["Masses"][:][order].astype(np.float64),
            Z=gas["MetalMassFractions"][:, -1][order].astype(np.float64),
            u={
                "FUV": gas["FUVSpecificEnergies"][:][order].astype(np.float64),
                "LW": gas["LWSpecificEnergies"][:][order].astype(np.float64),
            },
            F={
                "FUV": gas["FUVSpecificFluxes"][:][order].astype(np.float64),
                "LW": gas["LWSpecificFluxes"][:][order].astype(np.float64),
            },
            star_pos=star["Coordinates"][:].astype(np.float64),
            star_h=star["SmoothingLengths"][:].astype(np.float64),
            L={
                "FUV": star["FUVLuminosities"][:].astype(np.float64),
                "LW": star["LWLuminosities"][:].astype(np.float64),
            },
        )
    snap["c"] = C_LIGHT_CGS * snap["ut"] / snap["ul"]
    return snap


def snapshots(run: str) -> list:
    """Return the sorted snapshot files of a run, one per output time."""
    files, times = [], set()
    for path in sorted(glob.glob(os.path.join(run, "snap", "snapshot_*.hdf5"))):
        with h5py.File(path, "r") as f:
            t = float(np.asarray(f["/Header"].attrs["Time"]).flat[0])
        if not any(abs(t - s) <= 1e-9 * max(abs(t), 1e-300) for s in times):
            times.add(t)
            files.append(path)
    if not files:
        raise RuntimeError(f"No snapshots in {run}/snap")
    return files


def step_table(run: str) -> np.ndarray:
    """Return timesteps.txt rows: time, dt, gas updates, star updates."""
    rows = []
    with open(os.path.join(run, "timesteps.txt")) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            v = line.split()
            rows.append([float(v[1]), float(v[4]), int(v[7]), int(v[9])])
    return np.array(rows)


def c_hyp_margin(run: str) -> float:
    """Read GEARFeedback:ISRF_c_hyp_margin from used_parameters.yml."""
    with open(os.path.join(run, "used_parameters.yml")) as f:
        m = re.search(r"ISRF_c_hyp_margin:\s*([0-9.eE+-]+)", f.read())
    return float(m.group(1))


def kappa_rho(snap: dict, band: str) -> np.ndarray:
    """Return 1/lambda per particle, internal units."""
    kappa_cgs = SIGMA_D_CGS[band] * (snap["Z"] / GRACKLE_SOLAR_Z) / (MU_H * M_H_CGS)
    rho_cgs = snap["rho"] * snap["um"] / snap["ul"] ** 3
    return kappa_cgs * rho_cgs * snap["ul"]


def wendland_c2(r: np.ndarray, support: float) -> np.ndarray:
    """Return the 3D Wendland C2 kernel value for a support radius."""
    q = r / support
    w = 21.0 / (2.0 * np.pi * support**3) * (1.0 - q) ** 4 * (4.0 * q + 1.0)
    return np.where(q < 1.0, w, 0.0)


def periodic_dx(a: np.ndarray, b: np.ndarray, boxsize: float) -> np.ndarray:
    """Return minimum-image separations a - b."""
    d = a - b
    return d - boxsize * np.round(d / boxsize)


def injected_power(snap: dict, band: str) -> float:
    """Return P = sum_stars sum_j w_j L ext_j, internal units."""
    kappa_cgs = SIGMA_D_CGS[band] * (snap["Z"] / GRACKLE_SOLAR_Z) / (MU_H * M_H_CGS)
    sigma_cgs = 2.0 * GAMMA_3D * snap["h"] * snap["rho"] * snap["um"] / snap["ul"] ** 2
    ext = np.exp(-kappa_cgs * sigma_cgs)
    tree = cKDTree(snap["pos"], boxsize=snap["boxsize"])
    power = 0.0
    for s in range(len(snap["star_h"])):
        support = GAMMA_3D * snap["star_h"][s]
        idx = np.array(tree.query_ball_point(snap["star_pos"][s], support))
        r = np.linalg.norm(
            periodic_dx(snap["pos"][idx], snap["star_pos"][s], snap["boxsize"]), axis=1
        )
        mw = snap["mass"][idx] * wendland_c2(r, support)
        power += snap["L"][band][s] * np.sum(mw * ext[idx]) / np.sum(mw)
    return power


def identity_error(snap: dict, band: str) -> float:
    """Return the relative error of the amplitude identity."""
    lhs = np.sum(snap["mass"] * snap["u"][band] * kappa_rho(snap, band))
    rhs = injected_power(snap, band) / snap["c"]
    return float((lhs - rhs) / rhs)


def reduced_flux(snap: dict, band: str, run: str) -> np.ndarray:
    """Return f = |F| / (c_hyp u), c_hyp = min(C h / dt, c)."""
    dt = step_table(run)[-1, 1]
    c_hyp = np.minimum(c_hyp_margin(run) * snap["h"] / dt, snap["c"])
    u = snap["u"][band]
    Fmag = np.linalg.norm(snap["F"][band], axis=1)
    return np.where(u > 0.0, Fmag / (c_hyp * np.where(u > 0.0, u, 1.0)), np.nan)


def gate(name: str, value: float, bar: float) -> bool:
    """Print and return one gate; non-finite values fail."""
    ok = bool(np.isfinite(value) and np.isfinite(bar) and value <= bar)
    print(f"  {name}: {value:.4e} <= {bar:.4e} -> {'PASS' if ok else 'FAIL'}")
    return ok


def check_single_bin(run: str, n_gas: int, n_star: int) -> bool:
    """Require every gas particle and every star active on every step."""
    table = step_table(run)
    ok = bool(np.all(table[:, 2] == n_gas) and np.all(table[:, 3] == n_star))
    print(
        f"  {run}: all {n_gas} gas and {n_star} stars active every step "
        f"-> {'PASS' if ok else 'FAIL'}"
    )
    return ok


def lattice_sums(points: np.ndarray, d: float, lam: float) -> tuple:
    """Continuum lattice sums of the two steady kernels, over their mean.

    Sources sit at (k + 1/2) d for integer k. Free streaming:
    e^(-r/lam) / (4 pi r^2); diffusive: 3 e^(-sqrt(3) r/lam) / (4 pi lam r).
    Both kernels integrate to lam, so the uniform mean is n_s lam, with
    n_s = 1/d^3. Images beyond r_cut = 3 lam are added as their uniform
    mean.

    Parameters
    ----------
    points : np.ndarray
        Positions, shape (n, 3).
    d : float
        Lattice spacing.
    lam : float
        Absorption length.

    Returns
    -------
    tuple
        (free streaming, diffusive) fields over the mean, shape (n,) each.
    """
    r_cut = 3.0 * lam
    n = int(np.ceil(r_cut / d)) + 1
    k = np.arange(-n, n + 1)
    grid = (np.stack(np.meshgrid(k, k, k, indexing="ij"), -1).reshape(-1, 3) + 0.5) * d
    local = np.mod(points, d)
    fs = np.zeros(len(points))
    dif = np.zeros(len(points))
    for i, p in enumerate(local):
        r = np.linalg.norm(grid - p, axis=1)
        r = r[r < r_cut]
        fs[i] = np.sum(np.exp(-r / lam) / (4.0 * np.pi * r**2))
        dif[i] = np.sum(3.0 * np.exp(-np.sqrt(3.0) * r / lam) / (4.0 * np.pi * lam * r))
    tail_fs = np.exp(-r_cut / lam)
    x = np.sqrt(3.0) * r_cut / lam
    tail_dif = (1.0 + x) * np.exp(-x)
    mean = lam / d**3
    return fs / mean + tail_fs, dif / mean + tail_dif


def main() -> None:
    """Run every gate and exit nonzero on any failure."""
    opt = parse_options()
    run = lambda name: os.path.join(opt.dir, name)
    ok = True

    print("Time-bin preconditions")
    for name in ("A", "B", "AB", "injection_A", "injection_B", "injection_AB"):
        s = load(snapshots(run(name))[-1])
        ok &= check_single_bin(run(name), len(s["mass"]), len(s["star_h"]))

    print("E: exact injection superposition (propagation off, last snapshot)")
    inj = {k: load(snapshots(run("injection_" + k))[-1]) for k in ("A", "B", "AB")}
    eps = float(np.finfo(np.float32).eps)
    for band in ("FUV", "LW"):
        uA, uB, uAB = (inj[k]["u"][band] for k in ("A", "B", "AB"))
        total = uA + uB
        both = (uA > 0.0) & (uB > 0.0)
        lit = (total > 0.0) | (uAB > 0.0)
        floor = 4.0 * eps * (uA.max() + uB.max())
        excess = np.abs(uAB[lit] - total[lit]) - 1e-5 * total[lit]
        print(
            f"  {band}: {int(lit.sum())} illuminated, {int(both.sum())} reached "
            f"by both stars, {int(np.sum(lit & ((total == 0.0) | (uAB == 0.0))))} "
            f"lit in only one of AB and A+B"
        )
        worst = float(np.max(excess) / floor) if excess.size else np.nan
        ok &= gate(
            f"{band} max excess over 1e-5 relative, in float32 kernel floors",
            worst,
            1.0,
        )
        ok &= gate(
            f"{band} particles reached by both stars",
            1.0 / max(both.sum(), 1e-300),
            1.0,
        )
        mass = inj["AB"]["mass"]
        l1 = float(np.sum(mass * np.abs(uAB - total)) / np.sum(mass * total))
        print(f"  {band}: mass-weighted |u_AB - u_A - u_B| / (u_A + u_B) = {l1:.3e}")

    print("G1: amplitude identity with propagation on")
    series = {}
    for name in ("A", "B", "AB", "lattice", "lattice_single"):
        files = snapshots(run(name))[-opt.window :]
        snaps = [load(f) for f in files]
        series[name] = dict(
            last=snaps[-1],
            e={b: np.array([identity_error(s, b) for s in snaps]) for b in SIGMA_D_CGS},
        )
        print(
            f"  {name}: t = {snaps[-1]['time']:.4e}, e(FUV) over window = "
            f"{np.array2string(series[name]['e']['FUV'], precision=4)}, "
            f"e(LW) = {np.array2string(series[name]['e']['LW'], precision=4)}"
        )
    for band in ("FUV", "LW"):
        eA, eB, eAB = (series[k]["e"][band][-1] for k in ("A", "B", "AB"))
        sigma = max(np.std(series["A"]["e"][band]), np.std(series["B"]["e"][band]))
        ok &= gate(f"{band} |e_A|", abs(eA), 0.10)
        ok &= gate(f"{band} |e_B|", abs(eB), 0.10)
        ok &= gate(
            f"{band} |e_AB - (e_A + e_B)/2|",
            abs(eAB - 0.5 * (eA + eB)),
            0.5 * abs(eA - eB) + 3.0 * sigma,
        )

    print("Reported: float32 underflow of F.F (closure forced isotropic)")
    for name in ("A", "B", "AB", "lattice", "lattice_single"):
        s = series[name]["last"]
        for band in ("FUV", "LW"):
            Fmag = np.linalg.norm(s["F"][band], axis=1)
            print(
                f"  {name} {band}: min |F| = {Fmag.min():.3e}, below "
                f"{F_UNDERFLOW:.2e}: {int(np.sum(Fmag < F_UNDERFLOW))}"
            )

    print("Reported: propagated field of AB against u_A + u_B (no bar)")
    sA, sB, sAB = (series[k]["last"] for k in ("A", "B", "AB"))
    for band in ("FUV", "LW"):
        total = sA["u"][band] + sB["u"][band]
        share = sA["u"][band] / np.maximum(total, np.finfo(np.float64).tiny)
        rel = np.abs(sAB["u"][band] - total) / np.maximum(
            total, np.finfo(np.float64).tiny
        )
        dominated = (share > 0.9) | (share < 0.1)
        for label, sel in (("one source > 90%", dominated), ("mixed", ~dominated)):
            if sel.any():
                print(
                    f"  {band} {label}: {int(sel.sum())} particles, median "
                    f"{np.median(rel[sel]):.4e}, p90 {np.percentile(rel[sel], 90):.4e}"
                )
        for name, s in (("A", sA), ("AB", sAB)):
            print(
                f"  {band} {name}: median f = "
                f"{np.nanmedian(reduced_flux(s, band, run(name))):.3f}"
            )

    print("L1, L2: source lattice")
    lat = series["lattice"]["last"]
    one = series["lattice_single"]["last"]
    n_star = len(lat["star_h"])
    n_side = int(round(n_star ** (1.0 / 3.0)))
    d = lat["boxsize"] / n_side
    tree = cKDTree(lat["star_pos"], boxsize=lat["boxsize"])
    r_near, _ = tree.query(lat["pos"])
    outside = r_near > GAMMA_3D * np.max(lat["star_h"])
    print(
        f"  {n_star} stars, spacing {d:.4e}, {int(outside.sum())} gas outside kernels"
    )
    rng = np.random.default_rng(1)
    sample = rng.choice(
        np.flatnonzero(outside), min(opt.sample, int(outside.sum())), replace=False
    )
    r_one = np.linalg.norm(
        periodic_dx(one["pos"], one["star_pos"][0], one["boxsize"]), axis=1
    )
    for band in ("FUV", "LW"):
        e_lat = series["lattice"]["e"][band]
        e_one = series["lattice_single"]["e"][band]
        ok &= gate(
            f"{band} |e_lattice|", abs(e_lat[-1]), abs(e_one[-1]) + 3.0 * np.std(e_one)
        )
        f_lat = reduced_flux(lat, band, run("lattice"))
        ok &= gate(
            f"{band} median f outside kernels",
            float(np.nanmedian(f_lat[outside])),
            F_ISOTROPIC_BAR,
        )

        lam = 1.0 / np.median(kappa_rho(lat, band))
        u_mean = np.sum(lat["mass"] * lat["u"][band]) / np.sum(lat["mass"])
        scatter = np.sqrt(np.mean((lat["u"][band][sample] / u_mean - 1.0) ** 2))
        fs, dif = lattice_sums(lat["pos"][sample], d, lam)
        print(
            f"  {band}: d/lambda = {d / lam:.3f}; rms(u/u_mean - 1) outside "
            f"kernels: simulation {scatter:.4e}, continuum free streaming "
            f"{np.sqrt(np.mean((fs / np.mean(fs) - 1.0) ** 2)):.4e}, continuum "
            f"diffusive {np.sqrt(np.mean((dif / np.mean(dif) - 1.0) ** 2)):.4e}"
        )
        lo, hi = np.percentile(r_near[outside], [5, 95])
        band_sel = (r_one > lo) & (r_one < hi)
        f_one = reduced_flux(one, band, run("lattice_single"))
        print(
            f"  {band}: one-star run, median f at the same distances "
            f"[{lo:.3e}, {hi:.3e}]: {np.nanmedian(f_one[band_sel]):.3f}"
        )

    print("OVERALL:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
