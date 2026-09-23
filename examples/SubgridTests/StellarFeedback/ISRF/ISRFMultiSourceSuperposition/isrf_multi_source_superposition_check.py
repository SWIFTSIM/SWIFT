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

Gates (each per band, PE and LW):

E   Exact injection superposition, propagation off. Each illuminated gas
    particle gets an allowance
    a_i = 1e-5 (u(A)_i + u(B)_i) + 24 eps (u0(A) + u0(B)), eps the float32
    machine epsilon and u0 the deposit a star makes at the centre of its
    own kernel. The gate is the mass-weighted residual over the
    mass-weighted allowance,
    sum_i m_i |u(AB)_i - u(A)_i - u(B)_i| / sum_i m_i a_i, and the bar
    is 1.
    Second term: SWIFT evaluates the Wendland C2 polynomial
    (1 - x)^4 (1 + 4x) as 4x^5 - 15x^4 + 20x^3 - 10x^2 + 1 by Horner in
    float32 and clamps the result at zero. The last addition cancels two
    numbers of order ten, so the absolute error reaches 12 eps of the
    central value, measured over x in [0, 1] with this build's own flags
    (11.5 eps without FMA contraction, 8.1 eps with it). Each star's
    kernel is evaluated twice, once in the pair run and once in its own
    run, so the residual carries up to 24 eps of each star's central
    deposit. The error is absolute, so at the support radius, where the
    deposit is 1e-7 of the central value, it exceeds the deposit itself.
    First term: the three runs do not hold the same gas state after four
    steps, because the photoelectric heating of the field they differ in
    feeds back on the temperature and so on the positions and smoothing
    lengths. This term is an allowance, not a bound; a failure of this
    gate means a superposition defect, or that the drift grew past 1e-5
    relative.
    The statistic sums over particles instead of taking their maximum:
    only a few dozen particles are illuminated, so a maximum is decided by
    whether one of them lands in the cancellation region, and three
    identical runs of this fixture put the largest single-particle
    residual at 0.17, 1.7e-3 and 1.7e-3 of its allowance. Gated against
    an earlier floor of 4 eps of the peak field, those same three runs
    read 1.52, -9.2e-6 and -9.2e-6 against a bar of 1, so the second of
    them failed and the other two passed. A sum is bounded by the summed
    allowance by construction. The maximum is printed as a diagnostic.
    Gated with it: the number of particles reached by both stars, against
    half the lens volume of the two kernels times the gas number density,
    so that the sums above measure superposition and not two disjoint
    deposits.
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
    Every f entering that median must be finite, so an unlit particle
    fails the gate instead of being dropped from it.
L3  Sampling of the source lattice: the star kernel support spans at
    least 0.55 of the source spacing. At 0.299 (level 6 at n_side 8) the
    superposed field keeps its mean while its spatial variance exceeds
    the continuum lattice sum by an order of magnitude, so L1 and L2 are
    then measured on a field the scheme does not reproduce; at 0.598 (the
    shipped level 5 at n_side 8) they hold. The ratio at which that turns
    over has not been measured, and the bar is a guard placed just under
    the known-good point, not that onset. The kernels of a cubic lattice
    first overlap at (4 pi / 3) (H/d)^3 = 1, H/d = 0.620, which is where
    the ratio stops being a geometric statement about isolated sources;
    that is a derivation, not a measurement.

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
SIGMA_D_CGS = {"PE": 9e-22, "LW": 1.5e-21}
MU_H = 1.4
M_H_CGS = 1.6726219e-24
GRACKLE_SOLAR_Z = 0.01295
# Wendland C2 kernel support over smoothing length, 3D
KERNEL_NAME = "Wendland C2"
GAMMA_3D = 1.936492
C_LIGHT_CGS = 2.99792458e10
F_ISOTROPIC_BAR = 0.18
# Star kernel support over source spacing: known bad at 0.299, known good
# at the 0.598 of this example, turnover between them unmeasured, so the
# bar is a guard just under the known-good point.
KERNEL_SUPPORT_OVER_SPACING_BAR = 0.55
# Above this ratio the kernels cover the cubic lattice, the corner of a
# cell being sqrt(3)/2 of the spacing from the nearest source, and L2 is
# left with no gas outside them.
KERNEL_SUPPORT_OVER_SPACING_COVERED = 0.866
# Below this |F|, F.F underflows in float32 and the closure turns isotropic.
F_UNDERFLOW = 1.08e-19
# Absolute error of the float32 Horner evaluation of the Wendland C2
# polynomial, in machine epsilons of its central value, and the number of
# evaluations the injection superposition residual carries (see gate E).
KERNEL_HORNER_EPS = 12.0
KERNEL_EVALUATIONS_PER_STAR = 2
# Relative allowance for the gas states of the A, B and AB runs drifting
# apart over the four injection steps.
GAS_STATE_DRIFT = 1e-5
# Fraction of the continuum two-kernel lens count gate E requires, leaving
# room for the discreteness of the glass and for the two stars carrying
# slightly different smoothing lengths.
KERNEL_OVERLAP_FRACTION = 0.5


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


def read_extinction_path(snapshot_path: str) -> tuple:
    """Return the receiver-side extinction path the run used.

    Parameters
    ----------
    snapshot_path : str
        Any snapshot file of the run, used to locate its run directory.

    Returns
    -------
    tuple of (str, float)
        The mechanism name, and its path R in kernel support radii. R is
        meaningful only for ``constant_kernel_path``, where the column is
        ``R * kernel_gamma * h * rho``; it is NaN otherwise.
    """
    import yaml

    directory = os.path.dirname(os.path.dirname(snapshot_path))
    with open(os.path.join(directory, "used_parameters.yml")) as handle:
        used = yaml.safe_load(handle)["GEARFeedback"]
    if "ISRF_extinction_path" not in used:
        # A run archived before the key existed recorded no value at all, and
        # the path then in force was two kernel support radii.
        return "constant_kernel_path", 2.0
    name = used["ISRF_extinction_path"]
    if name == "constant_kernel_path":
        return name, float(used["ISRF_extinction_path_in_kernel_radii"])
    if name == "pair_separation":
        return name, float("nan")
    raise RuntimeError(
        f"GEARFeedback:ISRF_extinction_path {name!r} is not a length this "
        "check mirrors. Rerun the fixture with constant_kernel_path or "
        "pair_separation, or extend this check to that mechanism's own "
        "length."
    )


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
        scheme = f["/HydroScheme"].attrs
        name = scheme["Kernel function"]
        name = name.decode() if isinstance(name, bytes) else str(name)
        gamma = float(np.asarray(scheme["Kernel gamma"]).flat[0])
        if name != KERNEL_NAME or abs(gamma - GAMMA_3D) > 1e-5 * GAMMA_3D:
            raise RuntimeError(
                f"{path} was written by a {name} build of support over "
                f"smoothing length {gamma}. Every kernel weight and every "
                f"support radius here is the {KERNEL_NAME} one, so rebuild "
                f"with --with-kernel=wendland-C2."
            )
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
                "PE": gas["PESpecificEnergies"][:][order].astype(np.float64),
                "LW": gas["LWSpecificEnergies"][:][order].astype(np.float64),
            },
            F={
                "PE": gas["PESpecificFluxes"][:][order].astype(np.float64),
                "LW": gas["LWSpecificFluxes"][:][order].astype(np.float64),
            },
            star_pos=star["Coordinates"][:].astype(np.float64),
            star_h=star["SmoothingLengths"][:].astype(np.float64),
            L={
                "PE": star["PELuminosities"][:].astype(np.float64),
                "LW": star["LWLuminosities"][:].astype(np.float64),
            },
        )
    snap["c"] = C_LIGHT_CGS * snap["ut"] / snap["ul"]
    snap["ext_mechanism"], snap["ext_path_R"] = read_extinction_path(path)
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


def wendland_c2_polynomial(x: np.ndarray) -> np.ndarray:
    """Return the Wendland C2 polynomial (1 - x)^4 (1 + 4x), x = r / H."""
    return (1.0 - x) ** 4 * (1.0 + 4.0 * x)


def wendland_c2(r: np.ndarray, support: float) -> np.ndarray:
    """Return the 3D Wendland C2 kernel value for a support radius."""
    q = r / support
    w = 21.0 / (2.0 * np.pi * support**3) * wendland_c2_polynomial(q)
    return np.where(q < 1.0, w, 0.0)


def central_deposit(snap: dict, band: str) -> float:
    """Return the deposit a single star makes at the centre of its kernel.

    Parameters
    ----------
    snap : dict
        Snapshot of a run holding exactly one star, propagation off.
    band : str
        Band name.

    Returns
    -------
    float
        The largest u_i / W(x_i) over the particles inside 0.9 support
        radii, W the Wendland C2 polynomial and x_i = r_i / H. Outside
        that radius the polynomial itself is float32 noise, so those
        particles cannot report the central value. NaN if no particle sits
        inside it.
    """
    r = np.linalg.norm(
        periodic_dx(snap["pos"], snap["star_pos"][0], snap["boxsize"]), axis=1
    )
    x = r / (GAMMA_3D * snap["star_h"][0])
    inner = x < 0.9
    if not inner.any():
        return float("nan")
    return float(np.max(snap["u"][band][inner] / wendland_c2_polynomial(x[inner])))


def kernel_overlap_count(snap: dict) -> float:
    """Return the continuum count of gas inside both star kernels.

    Parameters
    ----------
    snap : dict
        Snapshot of the pair run, holding exactly two stars.

    Returns
    -------
    float
        The lens volume of the two kernel spheres times the mean gas
        number density, or 0 if the kernels do not overlap.
    """
    r1, r2 = GAMMA_3D * snap["star_h"][0], GAMMA_3D * snap["star_h"][1]
    d = float(
        np.linalg.norm(
            periodic_dx(snap["star_pos"][0], snap["star_pos"][1], snap["boxsize"])
        )
    )
    if d >= r1 + r2:
        return 0.0
    lens = (
        np.pi
        * (r1 + r2 - d) ** 2
        * (d**2 + 2.0 * d * (r1 + r2) - 3.0 * (r1 - r2) ** 2)
        / (12.0 * d)
    )
    return float(lens * len(snap["mass"]) / snap["boxsize"] ** 3)


def periodic_dx(a: np.ndarray, b: np.ndarray, boxsize: float) -> np.ndarray:
    """Return minimum-image separations a - b."""
    d = a - b
    return d - boxsize * np.round(d / boxsize)


def injected_power(snap: dict, band: str) -> float:
    """Return P = sum_stars sum_j w_j L ext_j, internal units."""
    kappa_cgs = SIGMA_D_CGS[band] * (snap["Z"] / GRACKLE_SOLAR_Z) / (MU_H * M_H_CGS)
    rho_cgs = snap["rho"] * snap["um"] / snap["ul"] ** 3
    pair_path = snap["ext_mechanism"] == "pair_separation"
    if not pair_path:
        ext = np.exp(
            -kappa_cgs
            * rho_cgs
            * snap["ext_path_R"]
            * GAMMA_3D
            * snap["h"]
            * snap["ul"]
        )
    tree = cKDTree(snap["pos"], boxsize=snap["boxsize"])
    power = 0.0
    for s in range(len(snap["star_h"])):
        support = GAMMA_3D * snap["star_h"][s]
        idx = np.array(tree.query_ball_point(snap["star_pos"][s], support))
        r = np.linalg.norm(
            periodic_dx(snap["pos"][idx], snap["star_pos"][s], snap["boxsize"]), axis=1
        )
        mw = snap["mass"][idx] * wendland_c2(r, support)
        # pair_separation shields each pair over its own separation, so the
        # factor cannot be hoisted out of the star loop.
        ext_j = (
            np.exp(-kappa_cgs[idx] * rho_cgs[idx] * r * snap["ul"])
            if pair_path
            else ext[idx]
        )
        power += snap["L"][band][s] * np.sum(mw * ext_j) / np.sum(mw)
    return power


def identity_error(snap: dict, band: str) -> float:
    """Return the relative error of the amplitude identity."""
    lhs = np.sum(snap["mass"] * snap["u"][band] * kappa_rho(snap, band))
    rhs = injected_power(snap, band) / snap["c"]
    return float((lhs - rhs) / rhs)


def reduced_flux(snap: dict, band: str, run: str) -> np.ndarray:
    """Return f = |F| / (c_hyp u), c_hyp = min(C h / dt, c).

    Parameters
    ----------
    snap : dict
        Snapshot fields.
    band : str
        Band name.
    run : str
        Run directory, for the last time step.

    Returns
    -------
    np.ndarray
        f per gas particle; an unlit particle gives NaN.
    """
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


def gate_at_least(name: str, value: float, bar: float) -> bool:
    """Print and return one lower-bound gate; non-finite values fail."""
    ok = bool(np.isfinite(value) and np.isfinite(bar) and value >= bar)
    print(f"  {name}: {value:.4e} >= {bar:.4e} -> {'PASS' if ok else 'FAIL'}")
    return ok


def median_gated(values: np.ndarray) -> float:
    """Return the median of a gated sample.

    Parameters
    ----------
    values : np.ndarray
        Sample to reduce.

    Returns
    -------
    float
        The median, or NaN if the sample is empty or holds a non-finite
        value, so that the gate fails instead of dropping it.
    """
    n_bad = int(np.sum(~np.isfinite(values)))
    if values.size == 0 or n_bad > 0:
        print(f"    {n_bad} of {values.size} values non-finite")
        return float("nan")
    return float(np.median(values))


def median_finite(values: np.ndarray) -> tuple:
    """Return the median over the finite values and the number of others.

    Parameters
    ----------
    values : np.ndarray
        Sample to reduce, for a printed quantity with no gate.

    Returns
    -------
    tuple
        (median over the finite values or NaN, count of non-finite ones).
    """
    finite = np.isfinite(values)
    n_bad = int(np.sum(~finite))
    if not finite.any():
        return float("nan"), n_bad
    return float(np.median(values[finite])), n_bad


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
    for band in ("PE", "LW"):
        uA, uB, uAB = (inj[k]["u"][band] for k in ("A", "B", "AB"))
        total = uA + uB
        both = (uA > 0.0) & (uB > 0.0)
        lit = (total > 0.0) | (uAB > 0.0)
        mass = inj["AB"]["mass"]
        print(
            f"  {band}: {int(lit.sum())} illuminated, {int(both.sum())} reached "
            f"by both stars, {int(np.sum(lit & ((total == 0.0) | (uAB == 0.0))))} "
            f"lit in only one of AB and A+B"
        )
        central = central_deposit(inj["A"], band) + central_deposit(inj["B"], band)
        allowance = (
            GAS_STATE_DRIFT * total[lit]
            + KERNEL_EVALUATIONS_PER_STAR * KERNEL_HORNER_EPS * eps * central
        )
        residual = np.abs(uAB[lit] - total[lit])
        budget = float(np.sum(mass[lit] * allowance))
        ok &= gate(
            f"{band} mass-weighted residual over its allowance",
            float(np.sum(mass[lit] * residual) / budget) if budget > 0.0 else np.nan,
            1.0,
        )
        ok &= gate_at_least(
            f"{band} particles reached by both stars",
            float(both.sum()),
            KERNEL_OVERLAP_FRACTION * kernel_overlap_count(inj["AB"]),
        )
        l1 = float(np.sum(mass * np.abs(uAB - total)) / np.sum(mass * total))
        worst = float(np.max(residual / allowance)) if residual.size else np.nan
        print(
            f"  {band}: mass-weighted |u_AB - u_A - u_B| / (u_A + u_B) = {l1:.3e}, "
            f"largest single-particle residual over its allowance = {worst:.3e}"
        )

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
            f"  {name}: t = {snaps[-1]['time']:.4e}, e(PE) over window = "
            f"{np.array2string(series[name]['e']['PE'], precision=4)}, "
            f"e(LW) = {np.array2string(series[name]['e']['LW'], precision=4)}"
        )
    for band in ("PE", "LW"):
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
        for band in ("PE", "LW"):
            Fmag = np.linalg.norm(s["F"][band], axis=1)
            print(
                f"  {name} {band}: min |F| = {Fmag.min():.3e}, below "
                f"{F_UNDERFLOW:.2e}: {int(np.sum(Fmag < F_UNDERFLOW))}"
            )

    print("Reported: propagated field of AB against u_A + u_B (no bar)")
    sA, sB, sAB = (series[k]["last"] for k in ("A", "B", "AB"))
    for band in ("PE", "LW"):
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
        for name, snap in (("A", sA), ("AB", sAB)):
            med, n_unlit = median_finite(reduced_flux(snap, band, run(name)))
            print(f"  {band} {name}: median f = {med:.3f}, {n_unlit} unlit")

    print("L1, L2: source lattice")
    lat = series["lattice"]["last"]
    one = series["lattice_single"]["last"]
    n_star = len(lat["star_h"])
    n_side = int(round(n_star ** (1.0 / 3.0)))
    d = lat["boxsize"] / n_side
    tree = cKDTree(lat["star_pos"], boxsize=lat["boxsize"])
    r_near, _ = tree.query(lat["pos"])
    outside = r_near > GAMMA_3D * np.max(lat["star_h"])
    n_outside = float(outside.sum())
    print(f"  {n_star} stars, spacing {d:.4e}, {int(n_outside)} gas outside kernels")
    support_over_spacing = float(GAMMA_3D * np.median(lat["star_h"]) / d)
    ok_sampling = gate_at_least(
        "L3 star kernel support over source spacing",
        support_over_spacing,
        KERNEL_SUPPORT_OVER_SPACING_BAR,
    )
    if not ok_sampling:
        print(
            f"    The star injection kernel spans {support_over_spacing:.3f} of "
            f"the {d:.4e} source spacing, under the "
            f"{KERNEL_SUPPORT_OVER_SPACING_BAR:.2f} this gate requires, a "
            "guard under the 0.598 this example is measured good at. At "
            "0.299 the superposed field keeps its mean to several "
            "significant figures while its spatial variance exceeds the "
            "continuum lattice sum by an order of magnitude; the ratio at "
            "which that turns over has not been measured. The cause is the "
            "kernel-to-spacing ratio itself and no fix exists, so the "
            "lattice gates below may measure a field the scheme does not "
            "reproduce. Run the lattice "
            "with a star kernel that reaches the neighbouring sources: a "
            "coarser gas resolution at this n_side, or more sources per side "
            f"at this resolution. Keep the ratio under "
            f"{KERNEL_SUPPORT_OVER_SPACING_COVERED:.3f}, above which the "
            "kernels cover the lattice and L2 is left with no gas outside "
            "them."
        )
    ok &= ok_sampling
    ok_outside = gate_at_least("L2 gas outside every star's kernel", n_outside, 1.0)
    if not ok_outside:
        print(
            "    Every particle sits inside a star kernel, so L2 has nothing "
            f"to measure. The support spans {support_over_spacing:.3f} of the "
            f"source spacing, at or above the "
            f"{KERNEL_SUPPORT_OVER_SPACING_COVERED:.3f} that covers a cubic "
            "lattice. Lower n_side, or raise the resolution level, until it "
            "stays under that."
        )
        print("OVERALL: FAIL")
        sys.exit(1)
    ok &= ok_outside
    rng = np.random.default_rng(1)
    sample = rng.choice(
        np.flatnonzero(outside), min(opt.sample, int(outside.sum())), replace=False
    )
    r_one = np.linalg.norm(
        periodic_dx(one["pos"], one["star_pos"][0], one["boxsize"]), axis=1
    )
    for band in ("PE", "LW"):
        e_lat = series["lattice"]["e"][band]
        e_one = series["lattice_single"]["e"][band]
        ok &= gate(
            f"{band} |e_lattice|", abs(e_lat[-1]), abs(e_one[-1]) + 3.0 * np.std(e_one)
        )
        f_lat = reduced_flux(lat, band, run("lattice"))
        ok &= gate(
            f"{band} median f outside kernels",
            median_gated(f_lat[outside]),
            F_ISOTROPIC_BAR,
        )

        lam = 1.0 / np.median(kappa_rho(lat, band))
        u_sample = lat["u"][band][sample]
        u_mean = np.sum(lat["mass"] * lat["u"][band]) / np.sum(lat["mass"])
        fs, dif = lattice_sums(lat["pos"][sample], d, lam)
        # Scatter and its two continuum references are all taken about the
        # mean over the SAME sample, so the three are one statistic. The
        # sample sits outside every kernel and is therefore colder than the
        # box, and that offset is a separate quantity: folding it into the
        # scatter inflates it by the gas the mask removed, not by any
        # departure from the lattice sums.
        scatter = np.sqrt(np.mean((u_sample / np.mean(u_sample) - 1.0) ** 2))
        offset = np.mean(u_sample) / u_mean - 1.0
        print(
            f"  {band}: d/lambda = {d / lam:.3f}; rms(u/u_sample_mean - 1) "
            f"outside kernels: simulation {scatter:.4e}, continuum free "
            f"streaming {np.sqrt(np.mean((fs / np.mean(fs) - 1.0) ** 2)):.4e}, "
            f"continuum diffusive "
            f"{np.sqrt(np.mean((dif / np.mean(dif) - 1.0) ** 2)):.4e}; "
            f"outside-kernel mean over box mean - 1: {offset:+.4e}"
        )
        lo, hi = np.percentile(r_near[outside], [5, 95])
        band_sel = (r_one > lo) & (r_one < hi)
        f_one = reduced_flux(one, band, run("lattice_single"))
        med_one, n_unlit = median_finite(f_one[band_sel])
        print(
            f"  {band}: one-star run, median f at the same distances "
            f"[{lo:.3e}, {hi:.3e}]: {med_one:.3f}, {n_unlit} unlit"
        )

    print("OVERALL:", "PASS" if ok else "FAIL")
    sys.exit(0 if ok else 1)


if __name__ == "__main__":
    main()
