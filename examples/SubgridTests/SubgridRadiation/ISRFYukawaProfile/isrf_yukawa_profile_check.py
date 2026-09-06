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
Tier 1: check that the Yukawa screened-diffusion propagation is a correct
integrator of its own intended equation. Bins gas particles by radius from
the central star at a converged snapshot, computes u(r)*r, and checks its
semi-log slope against the analytic Green's function's -1/lambda (Eq.
fuv-yukawa-green, theory/GEAR/Radiation/02_fuv_isrf.tex). lambda is computed
here from the run's own kappa_eff/Z/rho with an independent formula (not by
calling the C code under test), so this is a genuine check, not a
tautology.

Does not validate whether the Yukawa equation itself matches real
interstellar radiation transport -- that is a separate, physics-level
question (Tier 2), not a numerics one.
"""

import argparse
import glob
import os
import re
import sys

import h5py
import numpy as np

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# src/feedback/GEAR/radiation.h
SIGMA_D_FUV_CGS = 9e-22
SIGMA_D_LW_CGS = 1.5e-21
MU_H = 1.4
M_H_CGS = 1.6726219e-24
# Grackle's SolarMetalFractionByMass default, src/feedback/GEAR/radiation.h
GRACKLE_SOLAR_Z = 0.01295


def parse_options():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "-s",
        "--snapshot",
        default="snap/snapshot_*.hdf5",
        help="Glob pattern for snapshots to consider (default: %(default)s)",
    )
    parser.add_argument("--n-bins", type=int, default=25, help="Number of radial bins.")
    parser.add_argument(
        "--tol",
        type=float,
        default=0.3,
        help="Max allowed relative error on the fitted slope (default: %(default)s)",
    )
    parser.add_argument(
        "--output",
        default="isrf_yukawa_profile_check.png",
        help="Output plot filename.",
    )
    parser.add_argument(
        "--log",
        default="output.log",
        help="Run log to read the kernel-second-moment lambda correction "
        "from (default: %(default)s).",
    )
    return parser.parse_args()


# Matches feedback_props_print's own "LW/FUV Yukawa lambda correction
# (measured/analytic)        = <value>" line.
YUKAWA_LAMBDA_CORRECTION_RE = re.compile(
    r"LW/FUV Yukawa lambda correction.*=\s*([0-9.eE+-]+)"
)


def parse_yukawa_lambda_correction(log_path):
    """Read the kernel-second-moment lambda correction factor
    (lambda_measured/lambda_analytic for the naive decay-timescale
    correspondence) that feedback_props_print prints at start-up, once per
    run, for this build's own compiled-in kernel/eta_neighbours/hydro
    dimensionality. Never hardcode this number here: like W_min, it is
    kernel/eta-specific and must be read from the run's own log (see
    radiation_compute_yukawa_kernel_second_moment)."""
    if not os.path.exists(log_path):
        raise RuntimeError(
            f"Log file {log_path!r} not found -- cannot read the "
            "kernel-second-moment lambda correction factor."
        )
    with open(log_path, "r") as f:
        for line in f:
            m = YUKAWA_LAMBDA_CORRECTION_RE.search(line)
            if m:
                return float(m.group(1))
    raise RuntimeError(
        f"Could not find the 'LW/FUV Yukawa lambda correction' line in "
        f"{log_path!r} -- was LW_FUV_propagation actually enabled for "
        "this run?"
    )


def load_snapshot(path):
    with h5py.File(path, "r") as f:
        header = f["/Header"]
        time = float(np.asarray(header.attrs["Time"]).flat[0])
        boxsize = np.asarray(header.attrs["BoxSize"], dtype=float).flatten()[0]
        units = f["/Units"]
        unit_length_cgs = float(
            np.asarray(units.attrs["Unit length in cgs (U_L)"]).flat[0]
        )
        unit_mass_cgs = float(np.asarray(units.attrs["Unit mass in cgs (U_M)"]).flat[0])

        gas = f["/PartType0"]
        pos = gas["Coordinates"][:, :]
        rho = gas["Densities"][:]
        u_fuv = gas["FUVSpecificEnergies"][:]
        u_lw = gas["LWSpecificEnergies"][:]
        Z = gas["MetalMassFractions"][:, -1]

        star = f["/PartType4"]
        star_pos = star["Coordinates"][0, :]

    return dict(
        time=time,
        boxsize=boxsize,
        unit_length_cgs=unit_length_cgs,
        unit_mass_cgs=unit_mass_cgs,
        pos=pos,
        rho=rho,
        u_fuv=u_fuv,
        u_lw=u_lw,
        Z=Z,
        star_pos=star_pos,
    )


def radial_distance(pos, star_pos, boxsize):
    """Minimum-image radial distance from the star, for a periodic box."""
    dx = pos - star_pos
    dx -= boxsize * np.round(dx / boxsize)
    return np.sqrt(np.sum(dx**2, axis=1))


def analytic_lambda_cgs(Z, rho_internal, unit_length_cgs, unit_mass_cgs, sigma_d_cgs):
    """kappa_eff = sigma_d * (Z/Z_grackle_sun) / (mu_H * m_H); lambda =
    1/(kappa_eff*rho), independent of radiation_isrf.c's own formula."""
    rho_cgs = rho_internal * unit_mass_cgs / unit_length_cgs**3
    D_relative = Z / GRACKLE_SOLAR_Z
    kappa_eff_cgs = sigma_d_cgs * D_relative / (MU_H * M_H_CGS)
    return 1.0 / (kappa_eff_cgs * rho_cgs)


def fit_slope(r, u, r_min, r_max):
    """Semi-log fit of u(r)*r vs r over [r_min, r_max]; returns (slope, lambda)."""
    mask = (r > r_min) & (r < r_max) & (u > 0)
    if mask.sum() < 3:
        return None, None
    log_ur = np.log(u[mask] * r[mask])
    slope, intercept = np.polyfit(r[mask], log_ur, 1)
    return slope, (-1.0 / slope if slope < 0 else np.inf)


def main():
    opt = parse_options()
    files = sorted(glob.glob(opt.snapshot))
    if not files:
        raise RuntimeError(f"No snapshots match {opt.snapshot}")

    # Use the last snapshot: the run is designed to reach a converged
    # steady state well before it ends (see README).
    snap = load_snapshot(files[-1])

    r = radial_distance(snap["pos"], snap["star_pos"], snap["boxsize"])
    order = np.argsort(r)
    r, u_fuv, u_lw = r[order], snap["u_fuv"][order], snap["u_lw"][order]
    Z_mean = float(np.mean(snap["Z"]))
    rho_mean = float(np.mean(snap["rho"]))

    lambda_fuv_cgs = analytic_lambda_cgs(
        Z_mean,
        rho_mean,
        snap["unit_length_cgs"],
        snap["unit_mass_cgs"],
        SIGMA_D_FUV_CGS,
    )
    lambda_lw_cgs = analytic_lambda_cgs(
        Z_mean,
        rho_mean,
        snap["unit_length_cgs"],
        snap["unit_mass_cgs"],
        SIGMA_D_LW_CGS,
    )
    lambda_fuv = lambda_fuv_cgs / snap["unit_length_cgs"]
    lambda_lw = lambda_lw_cgs / snap["unit_length_cgs"]

    # radiation_snapshot_part_propagation now bakes the kernel-second-moment
    # correction into kappa_FUV/kappa_LW itself (see
    # radiation_compute_yukawa_kernel_second_moment's own doxygen and
    # 02_fuv_isrf.tex's fuv-discrete-correspondence), so the simulation's
    # realized screening length should already match the raw physical
    # lambda above. Read and print it only as an informational diagnostic;
    # it must not enter the pass/fail comparison.
    lambda_correction = parse_yukawa_lambda_correction(opt.log)

    # Bin edges: exclude the star's own kernel (r too small, where discrete
    # injection geometry dominates) and the outer quarter of the half-box
    # (periodic images start to bias the minimum-image distance there).
    half_box = snap["boxsize"] / 2.0
    r_min = 2.0 * (snap["boxsize"] / snap["pos"].shape[0] ** (1.0 / 3.0))
    r_max = 0.7 * half_box

    edges = np.linspace(r_min, r_max, opt.n_bins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    u_fuv_binned = np.full(opt.n_bins, np.nan)
    u_lw_binned = np.full(opt.n_bins, np.nan)
    for i in range(opt.n_bins):
        sel = (r >= edges[i]) & (r < edges[i + 1])
        if sel.sum() > 0:
            u_fuv_binned[i] = np.mean(u_fuv[sel])
            u_lw_binned[i] = np.mean(u_lw[sel])

    valid = ~np.isnan(u_fuv_binned) & ~np.isnan(u_lw_binned)
    slope_fuv, lambda_fuv_fit = fit_slope(
        centres[valid], u_fuv_binned[valid], r_min, r_max
    )
    slope_lw, lambda_lw_fit = fit_slope(
        centres[valid], u_lw_binned[valid], r_min, r_max
    )

    def report(band, lambda_fit, lambda_analytic):
        if lambda_fit is None:
            print(f"{band}: could not fit a slope (too few valid bins).")
            return False
        rel_err = abs(lambda_fit - lambda_analytic) / lambda_analytic
        status = "PASS" if rel_err < opt.tol else "FAIL"
        print(
            f"{band}: lambda_measured={lambda_fit:.4e}, "
            f"lambda_analytic={lambda_analytic:.4e}, "
            f"rel_err={rel_err:.3f} -> {status}"
        )
        return rel_err < opt.tol

    print(f"Snapshot: {files[-1]} (t={snap['time']:.4e})")
    print(f"Mean Z={Z_mean:.4e}, mean rho={rho_mean:.4e} (internal units)")
    print(f"Fit radial range: [{r_min:.4e}, {r_max:.4e}] (internal units)")
    print(
        f"Yukawa lambda correction (from {opt.log}): {lambda_correction:.4f} "
        f"(informational only: baked into kappa_FUV/kappa_LW internally, "
        f"not applied here)"
    )
    ok_fuv = report("FUV", lambda_fuv_fit, lambda_fuv)
    ok_lw = report("LW", lambda_lw_fit, lambda_lw)

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.semilogy(
        centres[valid], u_fuv_binned[valid] * centres[valid], "o-", label="FUV: u(r)*r"
    )
    ax.semilogy(
        centres[valid], u_lw_binned[valid] * centres[valid], "s-", label="LW: u(r)*r"
    )
    r_line = np.linspace(r_min, r_max, 100)
    ax.semilogy(
        r_line,
        (u_fuv_binned[valid][0] * centres[valid][0])
        * np.exp(-(r_line - centres[valid][0]) / lambda_fuv),
        "--",
        color="C0",
        label=r"FUV analytic $\propto e^{-r/\lambda}$",
    )
    ax.semilogy(
        r_line,
        (u_lw_binned[valid][0] * centres[valid][0])
        * np.exp(-(r_line - centres[valid][0]) / lambda_lw),
        "--",
        color="C1",
        label=r"LW analytic $\propto e^{-r/\lambda}$",
    )
    ax.set_xlabel("r (internal length units)")
    ax.set_ylabel(r"$u(r) \cdot r$")
    ax.legend()
    fig.tight_layout()
    fig.savefig(opt.output, dpi=150)
    print(f"Plot saved to {opt.output}")

    if not (ok_fuv and ok_lw):
        sys.exit(1)


if __name__ == "__main__":
    main()
