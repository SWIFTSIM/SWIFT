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
Tier 1: check that the hyperbolic (Cattaneo/P1-relaxation) LW/FUV
propagation scheme is a correct integrator of its OWN discretized
equations.

The governing PDE's steady state, in the continuum, is Yukawa-shaped
(`u(r) = S/(4*pi*D) * exp(-r/lambda)/r`), but this script does not compare
against that continuum profile. The discrete SPH estimator's actual fixed
point departs from the continuum profile whenever `h` is not much smaller
than `lambda` -- true at this project's own production resolution -- and a
continuum-target comparison reports a spurious failure there (33%-249%
error, `theory/GEAR/Radiation/verify_design_b_discrete_steady_state_
production_corners.py` and `.claude/dev/logs/
2026-09-08_1107_design-b-discrete-steady-state-vs-production-corners.md`)
even when the C code is solving its own discretized equations exactly.

Instead, this script builds the DISCRETE steady-state prediction directly:
it takes the run's own actual particle positions, smoothing lengths,
densities and masses from the snapshot, assembles the exact pairwise
`div(F)`/`grad(u)` estimators of `design-lw-fuv-design-b.md` Sec. 2.2, and
iterates Sec. 4.5's staggered exact-relaxation update to ITS OWN fixed
point (matrix-free Jacobi-style iteration; no linear solve needed, and the
fixed point is provably independent of `dt`/`c_hyp`, so an arbitrary,
numerically-convenient `c_hyp`/`dt` pair is used for the iteration only).
The measured simulation profile is then compared against this real-glass
discrete prediction, not the continuum Yukawa profile.

Tolerance: the discrete-solve methodology above was validated against four
actual production-resolution runs (`h/lambda` = 2.82 to 19.93, the
project's own thin-screening through deep-clump corners) in the
2026-09-08 log cited above, and matched the real simulation's fitted
screening length to 0.9%-6.2% at every corner, with the residual growing
(but not blowing up) with `h/lambda`. The default `--tol=0.15` (15%) is
about 2.4x the worst residual measured in that validation (6.2%), a
margin for run-to-run glass-disorder noise and metallicity/density
scatter this example's own corner was not part of. It is NOT a tolerance
chosen to make a marginal case pass; a script failure at this tolerance
should be treated as a real regression, not adjusted away.

Does not validate whether the governing PDE itself matches real
interstellar radiation transport -- that is a separate, physics-level
question (Tier 2), not a numerics one.

**The plotted "discrete-solve prediction" curve's ABSOLUTE amplitude is
not physically calibrated and must never be compared visually against the
simulation curve's absolute height.** `discrete_steady_state_lambda`
solves the fixed point for an arbitrary unit-total source
(`S /= np.sum(S)`, not the star's real `L_FUV`/`L_LW`, receiver-side
extinction, or the `3*c_hyp/c` rescale of `design-lw-fuv-design-b.md`
Sec 1.2), because only the SLOPE (the fitted `lambda`) is used for the
pass/fail gate below -- this is `design-lw-fuv-design-b.md` Sec 6.1's own
documented "amplitude criterion... still not implemented" gap, confirmed
empirically 2026-09-08 (see `.claude/dev/logs/
2026-09-08_1408_design-b-amplitude-units-investigation.md`): the two
curves in the saved plot can legitimately differ by 8-9 orders of
magnitude in absolute
height with no bug present, purely because of this normalization choice.
That investigation independently verified, from the star's actual
`L_FUV`/`L_LW`, the run's actual kappa/lambda, and the actual `3*c_hyp/c`
rescale, that the SIMULATION's own absolute amplitude is correct to
better than 0.01% (the exact identity `sum_i(m_i u_i) = (3/c) *
sum_over_star_kernel_neighbours(weight_j * lambda_j * L * extinction_j)`,
`c_hyp` cancelling algebraically per §1.2) -- the absolute-amplitude gap
visible in this script's own plot is a property of this CHECK's chosen
normalization, not of the code being checked.
"""

import argparse
import glob
import sys

import h5py
import matplotlib
import numpy as np
from scipy.spatial import cKDTree

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# src/feedback/GEAR/radiation.h
SIGMA_D_FUV_CGS = 9e-22
SIGMA_D_LW_CGS = 1.5e-21
MU_H = 1.4
M_H_CGS = 1.6726219e-24
# Grackle's SolarMetalFractionByMass default, src/feedback/GEAR/radiation.h
GRACKLE_SOLAR_Z = 0.01295
# Wendland C2, 3D: src/kernel_hydro.h kernel_gamma (this project's own
# --with-kernel=wendland-C2 build); needed to reconstruct the same kernel
# support radius (H = GAMMA_3D * h) the C code itself uses.
GAMMA_3D = 1.936492


def parse_options():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
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
        default=0.15,
        help="Max allowed relative error of the measured screening length "
        "against the discrete-solve prediction (default: %(default)s; "
        "see this script's own docstring for the justification).",
    )
    parser.add_argument(
        "--iter-tol",
        type=float,
        default=1e-6,
        help="Relative-change convergence tolerance for the discrete fixed-"
        "point iteration (default: %(default)s).",
    )
    parser.add_argument(
        "--max-iter",
        type=int,
        default=4000,
        help="Maximum iterations for the discrete fixed-point solve "
        "(default: %(default)s; converges in a few hundred iterations at "
        "every corner tested so far).",
    )
    parser.add_argument(
        "--output",
        default="isrf_hyperbolic_propagation_check.png",
        help="Output plot filename.",
    )
    return parser.parse_args()


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
        rho = gas["Densities"][:].astype(np.float64)
        h = gas["SmoothingLengths"][:].astype(np.float64)
        mass = gas["Masses"][:].astype(np.float64)
        u_fuv = gas["FUVSpecificEnergies"][:].astype(np.float64)
        u_lw = gas["LWSpecificEnergies"][:].astype(np.float64)
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
        h=h,
        mass=mass,
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
    1/(kappa_eff*rho): the physical dust-screening length, independent of
    this project's own kernel/eta_neighbours choice -- see
    radiation_get_part_linear_absorption_rate in radiation_isrf.c."""
    rho_cgs = rho_internal * unit_mass_cgs / unit_length_cgs**3
    D_relative = Z / GRACKLE_SOLAR_Z
    kappa_eff_cgs = sigma_d_cgs * D_relative / (MU_H * M_H_CGS)
    return 1.0 / (kappa_eff_cgs * rho_cgs)


def fit_slope(r, u, r_min, r_max):
    """Semi-log fit of u(r)*r vs r over [r_min, r_max]; returns (slope,
    lambda). Identical methodology for the real simulation output and the
    discrete-solve prediction, so the two are compared on equal footing."""
    mask = (r > r_min) & (r < r_max) & (u > 0)
    if mask.sum() < 3:
        return None, None
    log_ur = np.log(u[mask] * r[mask])
    slope, intercept = np.polyfit(r[mask], log_ur, 1)
    return slope, (-1.0 / slope if slope < 0 else np.inf)


# -----------------------------------------------------------------------------
# The discrete steady-state solve: real particle positions/h/rho/mass, exact
# Sec 2.2 pairwise div(F)/grad(u) estimators, Sec 4.5 staggered exact-
# relaxation iteration run to its own fixed point. Ported from the
# investigation script `theory/GEAR/Radiation/
# verify_design_b_discrete_steady_state_production_corners.py` (Part 5),
# validated there against four real production-resolution runs.
# -----------------------------------------------------------------------------


def wc2_3d_dwdr(r, H):
    """dW/dr of the 3D Wendland C2 kernel, support radius H (per-element
    array, one per pair-side, matching each particle's own h)."""
    H = np.broadcast_to(np.asarray(H, dtype=float), np.shape(r))
    q = r / H
    inside = q < 1.0
    out = np.zeros_like(r, dtype=float)
    qi, Hi = q[inside], H[inside]
    norm_i = 21.0 / (2.0 * np.pi * Hi**3)
    out[inside] = (
        norm_i * (-4.0 * (1.0 - qi) ** 3 * (4.0 * qi + 1.0) + 4.0 * (1.0 - qi) ** 4) / Hi
    )
    return out


def wc2_3d_w(r, H):
    """W of the 3D Wendland C2 kernel, support radius H."""
    H = np.broadcast_to(np.asarray(H, dtype=float), np.shape(r))
    q = r / H
    inside = q < 1.0
    out = np.zeros_like(r, dtype=float)
    qi, Hi = q[inside], H[inside]
    out[inside] = 21.0 / (2.0 * np.pi * Hi**3) * (1.0 - qi) ** 4 * (4.0 * qi + 1.0)
    return out


def build_pairs(pos, h, boxsize):
    """Periodic neighbour pairs (i<j) within either particle's own kernel
    support radius, plus dx, r and each side's own dW/dr."""
    H = GAMMA_3D * h
    tree = cKDTree(pos, boxsize=boxsize)
    pairs = np.array(sorted(tree.query_pairs(r=float(H.max()))))
    ii, jj = pairs[:, 0], pairs[:, 1]
    dx = pos[ii] - pos[jj]
    dx -= boxsize * np.round(dx / boxsize)
    r = np.linalg.norm(dx, axis=1)
    keep = (r > 0) & ((r < H[ii]) | (r < H[jj]))
    ii, jj, dx, r = ii[keep], jj[keep], dx[keep], r[keep]
    wi_dr = wc2_3d_dwdr(r, H[ii])
    wj_dr = wc2_3d_dwdr(r, H[jj])
    return ii, jj, dx, r, wi_dr, wj_dr


def grad_u(u, ii, jj, dx, r, wi_dr, wj_dr, rho, mass):
    """diffmode==0: grad(rho*u)/rho^2, design-lw-fuv-design-b.md Sec 2.2."""
    d_ij = rho[ii] * u[ii] - rho[jj] * u[jj]
    rinv = 1.0 / r
    coef_i = -mass[jj] * d_ij * wi_dr * rinv / rho[ii] ** 2
    coef_j = -mass[ii] * d_ij * wj_dr * rinv / rho[jj] ** 2
    g = np.zeros((3, len(u)))
    for k in range(3):
        np.add.at(g[k], ii, coef_i * dx[:, k])
        np.add.at(g[k], jj, coef_j * dx[:, k])
    return g


def div_F(Fvec, ii, jj, dx, r, wi_dr, wj_dr, rho, mass):
    """diffmode==1: shared-coefficient divergence, design doc Sec 2.2."""
    rinv = 1.0 / r
    Fi_dot = np.einsum("nk,nk->n", Fvec[:, ii].T, dx)
    Fj_dot = np.einsum("nk,nk->n", Fvec[:, jj].T, dx)
    Phi = Fi_dot / rho[ii] * wi_dr * rinv + Fj_dot / rho[jj] * wj_dr * rinv
    d = np.zeros(len(rho))
    np.add.at(d, ii, mass[jj] * Phi)
    np.add.at(d, jj, -mass[ii] * Phi)
    return d


def phi_relaxation_factor(a):
    """Identical to radiation_relaxation_phi_factor() in radiation_isrf.c:
    phi = (1-exp(-a))/a, series-expanded below a ~ 1e-6."""
    return np.where(a < 1e-6, 1.0 - 0.5 * a, -np.expm1(-a) / a)


def discrete_steady_state_lambda(snap, lam, r_min, r_max, n_bins, n_iter, iter_tol):
    """Solve Sec 4.5's staggered exact-relaxation update to its own fixed
    point on the run's REAL particle positions/h/rho/mass, for one band's
    physical screening length `lam`, then fit lambda_eff with the same
    methodology as the real simulation output (fit_slope above). `c_hyp`
    and `dt` are set to an arbitrary, numerically convenient value: the
    fixed point of this iteration is independent of both (design doc Sec
    4.5, 'Asymptotic preserving'; verified in `verify_design_b_
    timestepping_stability.py` Part C.2, profiles at dt/tau = 0.004, 0.04,
    0.2 agree to 1e-8), so no attempt is made to match the real run's own
    per-step dt/c_hyp history -- only the fixed point matters here, not
    the path to it.
    """
    pos, h, rho, mass = snap["pos"], snap["h"], snap["rho"], snap["mass"]
    boxsize = snap["boxsize"]
    N = len(pos)
    ii, jj, dx, r, wi_dr, wj_dr = build_pairs(pos, h, boxsize)

    # Kernel-weighted deposit at the star's position: an idealised stand-in
    # for the star's real injection-neighbour weighting (radiation_iact.h),
    # approximated here with the run's own median smoothing length. This is
    # the dominant known idealisation left in this methodology (see the
    # 2026-09-08 validation log's "Open questions"); it does not prevent
    # the 0.9%-6.2% agreement this script's tolerance is built around.
    dxs = pos - snap["star_pos"]
    dxs -= boxsize * np.round(dxs / boxsize)
    rs = np.linalg.norm(dxs, axis=1)
    h_med = np.median(h)
    Hs = GAMMA_3D * h_med
    S = wc2_3d_w(rs, Hs)
    S /= np.sum(S)

    # Arbitrary numerically-convenient closure (see docstring above): the
    # fixed point does not depend on this choice.
    C_HYP_ARBITRARY = 0.1
    c_hyp = C_HYP_ARBITRARY * h_med
    a = C_HYP_ARBITRARY * h_med / lam
    e_, ph = np.exp(-a), phi_relaxation_factor(a)

    u = np.zeros(N)
    F = np.zeros((3, N))
    it = 0
    du = np.inf
    for it in range(n_iter):
        dF = div_F(F, ii, jj, dx, r, wi_dr, wj_dr, rho, mass)
        u_new = e_ * u + ph * (S - dF)
        g = grad_u(u_new, ii, jj, dx, r, wi_dr, wj_dr, rho, mass)
        F_new = e_ * F - ph * c_hyp**2 * g
        du = np.max(np.abs(u_new - u)) / max(np.max(np.abs(u_new)), 1e-300)
        u, F = u_new, F_new
        if du < iter_tol and it > 20:
            break
    converged = du < iter_tol

    r_star = np.linalg.norm(dxs, axis=1)
    order = np.argsort(r_star)
    r_sorted, u_sorted = r_star[order], u[order]
    edges = np.linspace(r_min, r_max, n_bins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    u_binned = np.full(n_bins, np.nan)
    for i in range(n_bins):
        sel = (r_sorted >= edges[i]) & (r_sorted < edges[i + 1])
        if sel.sum() > 0:
            u_binned[i] = np.mean(u_sorted[sel])
    valid = ~np.isnan(u_binned)
    _, lam_fit = fit_slope(centres[valid], u_binned[valid], r_min, r_max)
    return lam_fit, centres[valid], u_binned[valid], it + 1, converged


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
    h_mean = float(np.median(snap["h"]))

    lambda_fuv_cgs = analytic_lambda_cgs(
        Z_mean, rho_mean, snap["unit_length_cgs"], snap["unit_mass_cgs"], SIGMA_D_FUV_CGS
    )
    lambda_lw_cgs = analytic_lambda_cgs(
        Z_mean, rho_mean, snap["unit_length_cgs"], snap["unit_mass_cgs"], SIGMA_D_LW_CGS
    )
    lambda_fuv = lambda_fuv_cgs / snap["unit_length_cgs"]
    lambda_lw = lambda_lw_cgs / snap["unit_length_cgs"]

    # Bin edges: exclude the star's own kernel (r too small, where discrete
    # injection geometry dominates) and the outer quarter of the half-box
    # (periodic images start to bias the minimum-image distance there).
    # Identical range for the real simulation fit and the discrete-solve
    # fit, so the two are compared like-for-like.
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

    slope_fuv, lambda_fuv_measured = fit_slope(
        centres[valid], u_fuv_binned[valid], r_min, r_max
    )
    slope_lw, lambda_lw_measured = fit_slope(
        centres[valid], u_lw_binned[valid], r_min, r_max
    )

    print(f"Snapshot: {files[-1]} (t={snap['time']:.4e})")
    print(f"Mean Z={Z_mean:.4e}, mean rho={rho_mean:.4e} (internal units)")
    print(f"h/lambda: FUV={h_mean / lambda_fuv:.3f}, LW={h_mean / lambda_lw:.3f}")
    print(f"Fit radial range: [{r_min:.4e}, {r_max:.4e}] (internal units)")
    print()

    print("Solving each band's discrete steady state on the run's own real")
    print("particle positions/h/rho/mass (may take a few seconds)...")
    lambda_fuv_discrete, dc_fuv, ud_fuv, it_fuv, conv_fuv = discrete_steady_state_lambda(
        snap, lambda_fuv, r_min, r_max, opt.n_bins, opt.max_iter, opt.iter_tol
    )
    lambda_lw_discrete, dc_lw, ud_lw, it_lw, conv_lw = discrete_steady_state_lambda(
        snap, lambda_lw, r_min, r_max, opt.n_bins, opt.max_iter, opt.iter_tol
    )

    def report(band, lambda_measured, lambda_discrete, lambda_continuum, iters, converged):
        if lambda_measured is None:
            print(f"{band}: could not fit the simulation's own slope (too few valid bins).")
            return False
        if lambda_discrete is None or not np.isfinite(lambda_discrete):
            print(f"{band}: discrete-solve prediction could not be fit either.")
            return False
        if not converged:
            print(
                f"{band}: WARNING -- discrete fixed-point iteration did not "
                f"converge within {iters} iterations; result may be unreliable."
            )
        rel_err = abs(lambda_measured - lambda_discrete) / lambda_discrete
        status = "PASS" if rel_err < opt.tol else "FAIL"
        print(
            f"{band}: lambda_measured(sim)={lambda_measured:.4e}, "
            f"lambda_discrete(prediction)={lambda_discrete:.4e}, "
            f"lambda_continuum(analytic)={lambda_continuum:.4e}, "
            f"rel_err(sim vs. discrete)={rel_err:.3f} -> {status} "
            f"[discrete solve: {iters} iterations, "
            f"{'converged' if converged else 'NOT converged'}]"
        )
        return rel_err < opt.tol

    ok_fuv = report("FUV", lambda_fuv_measured, lambda_fuv_discrete, lambda_fuv, it_fuv, conv_fuv)
    ok_lw = report("LW", lambda_lw_measured, lambda_lw_discrete, lambda_lw, it_lw, conv_lw)

    # Informational only (not part of the pass/fail gate): a global
    # negative-value count and the largest single non-monotonic bin-to-bin
    # increase in the radially-binned profile. Small near-injection-kernel
    # ripples and sub-per-mille negative undershoots are ordinary SPH
    # discreteness noise at this particle count, not validated here as a
    # separate leg -- see the README for what this script's own pass/fail
    # criterion actually is.
    n_negative = int(np.sum(snap["u_fuv"] < 0) + np.sum(snap["u_lw"] < 0))
    worst_bump_fuv = float(np.max(np.diff(u_fuv_binned[valid]) / u_fuv_binned[valid][:-1]))
    worst_bump_lw = float(np.max(np.diff(u_lw_binned[valid]) / u_lw_binned[valid][:-1]))
    print(
        f"(informational) particles with u < 0: {n_negative}/{2 * len(snap['u_fuv'])}; "
        f"largest bin-to-bin fractional increase in u(r): "
        f"FUV={worst_bump_fuv:.3f}, LW={worst_bump_lw:.3f}"
    )

    fig, ax = plt.subplots(figsize=(6, 5))
    ax.semilogy(centres[valid], u_fuv_binned[valid] * centres[valid], "o-", label="FUV: sim")
    ax.semilogy(centres[valid], u_lw_binned[valid] * centres[valid], "s-", label="LW: sim")
    ax.semilogy(
        dc_fuv, ud_fuv * dc_fuv, "--", color="C0",
        label="FUV: discrete-solve prediction (slope only -- see note below)",
    )
    ax.semilogy(
        dc_lw, ud_lw * dc_lw, "--", color="C1",
        label="LW: discrete-solve prediction (slope only -- see note below)",
    )
    ax.set_xlabel("r (internal length units)")
    ax.set_ylabel(r"$u(r) \cdot r$")
    ax.legend(fontsize=7, loc="upper right")
    ax.text(
        0.02, 0.02,
        "Dashed curves: unit-normalized source -> compare SLOPE only.\n"
        "Absolute height is not physically calibrated (see script docstring).",
        transform=ax.transAxes, fontsize=6.5, va="bottom", ha="left",
        color="0.4",
    )
    fig.tight_layout()
    fig.savefig(opt.output, dpi=150)
    print(f"Plot saved to {opt.output}")
    print(
        "NOTE: the dashed 'discrete-solve prediction' curves use an arbitrary "
        "unit-total source normalization (see this script's own docstring); "
        "their absolute height is not meant to match the solid simulation "
        "curves and any such gap (however large) is not, by itself, evidence "
        "of a bug. Only the fitted lambda (reported above) is a validated "
        "comparison."
    )

    if not (ok_fuv and ok_lw):
        sys.exit(1)


if __name__ == "__main__":
    main()
