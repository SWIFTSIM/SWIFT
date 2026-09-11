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
Tier 1: check that the hyperbolic (Cattaneo-type, M1-closure) LW/FUV
propagation scheme is a correct integrator of its OWN discretized
equations.

The governing equations, the closure, the flux limiter and both analytic
steady-state limits are derived in `theory/GEAR/Radiation/02_fuv_isrf.tex`,
Sections `fuv-p1` (closure and limiter), `fuv-lambda-chyp` (the
RSOL-consistent constants) and `fuv-steady-limits` (the two limits, worked
through). This script's own prediction is the discrete image of those
equations; read that section first, the algebra here is not repeated.

WHAT IS COMPARED, AND WHY NOT A CLOSED-FORM PROFILE
---------------------------------------------------
The continuum steady state of the implemented scheme is not a single
closed-form profile. It has two branches, selected locally by the reduced
flux `f = |F|/(c_hyp*u)`:

  free streaming (f -> 1):   u(r) = L*exp(-r/lambda) / (4*pi*c*r**2)
  diffusion     (f -> 0):   u(r) = 3*L*exp(-sqrt(3)*r/lambda) / (4*pi*c*lambda*r)

(theory doc Eqs. `fuv-steady-freestream`/`fuv-steady-diffusion`; both are
exact, and `c_hyp` cancels out of both). A single point source in a purely
absorbing medium sits on the free-streaming branch: substituting the
diffusive profile back into `f` gives `f -> 1/sqrt(3)`, not 0, so that
branch is not self-consistent here. The diffusive branch is the one a
many-source diffuse background lives on, which is this module's production
use case but not this example's configuration.

On top of that, the DISCRETE fixed point departs from either continuum
profile whenever `h` is not much smaller than `lambda` -- true at this
project's own production resolution -- so a continuum-target comparison
reports a spurious failure (33%-249% error,
`theory/GEAR/Radiation/verify_design_b_discrete_steady_state_production_
corners.py`) even when the C code is solving its own discretized equations
exactly.

So, as before, this script builds the DISCRETE steady-state prediction
directly: it takes the run's own actual particle positions, smoothing
lengths, densities and masses from the snapshot, assembles the exact
pairwise operators the C code uses (`radiation_propagation_iact.h`'s
`div(F)` and anisotropic-M1 `grad(u)`, plus `radiation_isrf.c`'s flux
limiter), and iterates the staggered exact-relaxation update to ITS OWN
fixed point. The measured simulation profile is then compared against that
real-glass discrete prediction.

THE PASS/FAIL CRITERION CHANGED WITH THE CLOSURE, DELIBERATELY
--------------------------------------------------------------
The P1-era version of this script gated on the fitted screening length
alone: it fitted `log(u*r)` against `r` for both the simulation and the
prediction and compared the two slopes. Under the Yukawa steady state that
fit is exact and its slope IS `-1/lambda_eff`. Under M1 it is not: on the
free-streaming branch `u*r ~ exp(-r/lambda)/r`, so the fitted slope picks
up a geometric `-log(r)` contribution that is identical in the simulation
and in the prediction by construction (same positions, same bins). A real
closure error can therefore hide behind a shared `1/r**2`, and a
lambda-only gate is no longer a trustworthy regression gate.

This script therefore gates on two closure-agnostic quantities instead,
and reports the fitted lambda as information only:

  1. `u(r)`, bin by bin, simulation against discrete prediction, over the
     same radial bins -- assumes no functional form at all, and subsumes
     both shape and amplitude. This is the primary gate (`--tol`).
  2. The total field amplitude, against a first-principles identity that
     needs no fit, no profile and no closure (theory doc
     Eq. `fuv-steady-amplitude`): summing the steady-state zeroth moment
     over all particles kills the divergence term exactly, because the
     discrete `div(F)` operator is a mirrored credit/debit pair with
     `sum_i m_i*(div F)_i = 0` identically. What remains is

         sum_i m_i*u_i/lambda_i = (1/c) * sum_j weight_j * L * extinction_j

     (exact for a uniform `c_hyp`; `weight_j`/`extinction_j` are the
     injection weights and receiver-side extinction of the theory doc's
     Eq. `fuv-inject`). Gated with a wider tolerance (`--amp-tol`), since
     the run's `c_hyp` is per-particle rather than exactly uniform.

The reduced flux `f` realized in the fit range is printed, so a reader can
see which branch the run is actually on rather than assuming one.

Does not validate whether the governing equations themselves match real
interstellar radiation transport -- that is a separate, physics-level
question (Tier 2), not a numerics one.
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
# Speed of light, cgs; the c_hyp/c source rescale of radiation_isrf.c's
# radiation_end_density_propagation.
C_LIGHT_CGS = 2.99792458e10


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
        help="Max allowed MEDIAN relative error of the binned u(r) profile "
        "against the discrete-solve prediction (default: %(default)s).",
    )
    parser.add_argument(
        "--amp-tol",
        type=float,
        default=0.10,
        help="Max allowed relative error of the total-amplitude identity "
        "(default: %(default)s).",
    )
    parser.add_argument(
        "--iter-tol",
        type=float,
        default=1e-6,
        help="Relative-change convergence tolerance for the discrete fixed-"
        "point iteration, on u AND on F (default: %(default)s).",
    )
    parser.add_argument(
        "--max-iter",
        type=int,
        default=4000,
        help="Maximum iterations for the discrete fixed-point solve "
        "(default: %(default)s).",
    )
    parser.add_argument(
        "--relax",
        type=float,
        default=1.0,
        help="Under-relaxation factor for the fixed-point iteration "
        "(default: %(default)s, i.e. none). Lower it if the iteration "
        "limit-cycles on the flux limiter's kink instead of converging; "
        "do NOT loosen --iter-tol instead.",
    )
    parser.add_argument(
        "--check-chyp-invariance",
        action="store_true",
        help="Re-solve the FUV band at a 5x SMALLER arbitrary c_hyp and "
        "report the profile difference. The discrete fixed point is "
        "c_hyp-independent analytically (see discrete_steady_state's "
        "docstring); this re-measures that property rather than assuming it. "
        "Smaller, not larger: the fixed point does not depend on c_hyp, but "
        "the iteration's own stability does, through the same Courant number "
        "c_hyp*dt/h the real scheme obeys.",
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
        unit_time_cgs = float(np.asarray(units.attrs["Unit time in cgs (U_t)"]).flat[0])

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
        star_h = float(star["SmoothingLengths"][0])
        L_FUV = float(star["FUVLuminosities"][0])
        L_LW = float(star["LWLuminosities"][0])

    return dict(
        time=time,
        boxsize=boxsize,
        unit_length_cgs=unit_length_cgs,
        unit_mass_cgs=unit_mass_cgs,
        unit_time_cgs=unit_time_cgs,
        pos=pos,
        rho=rho,
        h=h,
        mass=mass,
        u_fuv=u_fuv,
        u_lw=u_lw,
        Z=Z,
        star_pos=star_pos,
        star_h=star_h,
        L_FUV=L_FUV,
        L_LW=L_LW,
    )


def radial_distance(pos, star_pos, boxsize):
    """Minimum-image radial distance from the star, for a periodic box."""
    dx = pos - star_pos
    dx -= boxsize * np.round(dx / boxsize)
    return np.sqrt(np.sum(dx**2, axis=1))


def kappa_eff_mass_opacity_cgs(Z, sigma_d_cgs):
    """Band dust mass opacity (area/mass, cgs); mirrors radiation_get_dust_mass_opacity
    in radiation_isrf.c, with local_dust_to_gas_ratio left at Grackle's own
    default (this example's params.yml leaves it unset), so D_relative reduces
    to Z/Z_grackle_sun."""
    D_relative = Z / GRACKLE_SOLAR_Z
    return sigma_d_cgs * D_relative / (MU_H * M_H_CGS)


def analytic_lambda_cgs(Z, rho_internal, unit_length_cgs, unit_mass_cgs, sigma_d_cgs):
    """lambda = 1/(kappa_eff*rho): the physical dust-screening length,
    independent of this project's own kernel/eta_neighbours choice -- see
    radiation_get_part_linear_absorption_rate in radiation_isrf.c."""
    rho_cgs = rho_internal * unit_mass_cgs / unit_length_cgs**3
    kappa_eff_cgs = kappa_eff_mass_opacity_cgs(Z, sigma_d_cgs)
    return 1.0 / (kappa_eff_cgs * rho_cgs)


def receiver_extinction_factor(
    Z, rho_internal, h_internal, unit_length_cgs, unit_mass_cgs, sigma_d_cgs
):
    """Receiver-side dust extinction exp(-kappa_eff*Sigma_gas); mirrors
    radiation_get_part_LW_FUV_extinction_factors, with the comoving column
    density Sigma_gas = 2*kernel_gamma*h*rho (this example is
    non-cosmological, so comoving equals physical here)."""
    Sigma_gas_cgs = (
        2.0 * GAMMA_3D * h_internal * rho_internal * unit_mass_cgs / unit_length_cgs**2
    )
    kappa_eff_cgs = kappa_eff_mass_opacity_cgs(Z, sigma_d_cgs)
    return np.exp(-kappa_eff_cgs * Sigma_gas_cgs)


def fit_slope(r, u, r_min, r_max):
    """Semi-log fit of u(r)*r vs r over [r_min, r_max]; returns (slope,
    lambda). Informational only under the M1 closure (see this module's
    docstring): the same fit is applied to the simulation and to the
    discrete prediction, so the two remain comparable with each other, but
    the fitted number is not a screening length on the free-streaming
    branch."""
    mask = (r > r_min) & (r < r_max) & (u > 0)
    if mask.sum() < 3:
        return None, None
    log_ur = np.log(u[mask] * r[mask])
    slope, intercept = np.polyfit(r[mask], log_ur, 1)
    return slope, (-1.0 / slope if slope < 0 else np.inf)


def radial_bin(r, u, edges):
    """Mean of u in each radial bin; NaN where a bin is empty."""
    n_bins = len(edges) - 1
    out = np.full(n_bins, np.nan)
    for i in range(n_bins):
        sel = (r >= edges[i]) & (r < edges[i + 1])
        if sel.sum() > 0:
            out[i] = np.mean(u[sel])
    return out


# -----------------------------------------------------------------------------
# The discrete steady-state solve: real particle positions/h/rho/mass, the
# exact pairwise operators of src/feedback/GEAR/radiation_propagation_iact.h
# and the flux limiter of radiation_isrf.c, iterated to their own fixed
# point. Theory: theory/GEAR/Radiation/02_fuv_isrf.tex, Sections `fuv-p1`,
# `fuv-lambda-chyp`, `fuv-operators`, `fuv-steady-limits`.
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
        norm_i
        * (-4.0 * (1.0 - qi) ** 3 * (4.0 * qi + 1.0) + 4.0 * (1.0 - qi) ** 4)
        / Hi
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


def m1_closure_tensor(u, F, c_hyp):
    """Per-particle M1 closure tensor D(f), shape (N, 3, 3); mirrors
    radiation_get_m1_closure_tensor_band in radiation_propagation_iact.h,
    zero-flux guard included (n = 0 where |F| = 0, f = 0 where u <= 0, and
    f clamped at 1 before the sqrt(4-3f^2) is formed)."""
    F2 = np.einsum("kn,kn->n", F, F)
    F_inv = np.where(F2 > 0.0, 1.0 / np.sqrt(np.where(F2 > 0.0, F2, 1.0)), 0.0)
    Fmag = F2 * F_inv
    n = F * F_inv  # (3, N), exactly zero where |F| = 0

    denom = c_hyp * u
    f = np.where(
        denom > 0.0, np.minimum(Fmag / np.where(denom > 0.0, denom, 1.0), 1.0), 0.0
    )

    chi = (3.0 + 4.0 * f * f) / (5.0 + 2.0 * np.sqrt(4.0 - 3.0 * f * f))
    iso = 0.5 * (1.0 - chi)
    aniso = 0.5 * (3.0 * chi - 1.0)

    D = aniso[:, None, None] * np.einsum("an,bn->nab", n, n)
    idx = np.arange(3)
    D[:, idx, idx] += iso[:, None]
    return D, f


def grad_P(u, D, ii, jj, dx, r, wi_dr, wj_dr, rho, mass):
    """The anisotropic M1 pressure-tensor divergence accumulator,
    `(div(D*rho*u))/rho`, in the `diffmode == 0` (own kernel derivative,
    no shared average, no grad-h factor) form of
    radiation_gradient_accumulate_band. Reduces to the older scalar
    grad(rho*u)/(3*rho^2) form when D = I/3."""
    rinv = 1.0 / r
    Ui = rho[ii] * u[ii]
    Uj = rho[jj] * u[jj]
    # (D . dx) for each pair side, contracted on the tensor's second index.
    Di_dx = np.einsum("nab,nb->na", D[ii], dx)
    Dj_dx = np.einsum("nab,nb->na", D[jj], dx)
    T = (Ui * rinv)[:, None] * Di_dx - (Uj * rinv)[:, None] * Dj_dx

    fac_i = -mass[jj] * wi_dr / rho[ii] ** 2
    fac_j = -mass[ii] * wj_dr / rho[jj] ** 2

    g = np.zeros((3, len(u)))
    for k in range(3):
        np.add.at(g[k], ii, fac_i * T[:, k])
        np.add.at(g[k], jj, fac_j * T[:, k])
    return g


def div_F(Fvec, ii, jj, dx, r, wi_dr, wj_dr, rho, mass):
    """The shared-coefficient (`diffmode == 1`) divergence of
    radiation_divergence_accumulate_band. Closure-independent, unchanged by
    the P1-to-M1 upgrade."""
    rinv = 1.0 / r
    Fi_dot = np.einsum("nk,nk->n", Fvec[:, ii].T, dx)
    Fj_dot = np.einsum("nk,nk->n", Fvec[:, jj].T, dx)
    Phi = Fi_dot / rho[ii] * wi_dr * rinv + Fj_dot / rho[jj] * wj_dr * rinv
    d = np.zeros(len(rho))
    np.add.at(d, ii, mass[jj] * Phi)
    np.add.at(d, jj, -mass[ii] * Phi)
    return d


def apply_flux_limiter(u, F, c_hyp):
    """`F <- F*min(1, c_hyp*u/|F|)` for u > 0, `F <- 0` for u <= 0; mirrors
    radiation_apply_flux_limiter_band in radiation_isrf.c."""
    F2 = np.einsum("kn,kn->n", F, F)
    Fmag = np.sqrt(F2)
    limiter = np.ones_like(u)
    live = (u > 0.0) & (F2 > 0.0)
    limiter[live] = np.minimum(1.0, c_hyp * u[live] / Fmag[live])
    limiter[u <= 0.0] = 0.0
    return F * limiter


def phi_relaxation_factor(a):
    """Identical to radiation_relaxation_phi_factor() in radiation_isrf.c:
    phi = (1-exp(-a))/a, series-expanded below a ~ 1e-6."""
    return np.where(a < 1e-6, 1.0 - 0.5 * a, -np.expm1(-a) / a)


def injection_source(snap, sigma_d_band_cgs, L_band):
    """The run's own real injection footprint and receiver-side extinction:
    per-particle specific injection rate `S_true` (energy/mass/time, internal
    units), mirroring radiation_iact.h's `weight = mj*wi*si_inv_weight`.
    Also returns the total injected power `sum_j weight_j*L*extinction_j`,
    which the amplitude identity of the theory doc's
    Eq. `fuv-steady-amplitude` needs."""
    pos, mass, rho, h, Z = snap["pos"], snap["mass"], snap["rho"], snap["h"], snap["Z"]
    boxsize = snap["boxsize"]
    dxs = pos - snap["star_pos"]
    dxs -= boxsize * np.round(dxs / boxsize)
    rs = np.linalg.norm(dxs, axis=1)
    Hs = GAMMA_3D * snap["star_h"]
    mass_weighted_kernel = mass * wc2_3d_w(rs, Hs)
    weight = mass_weighted_kernel / np.sum(mass_weighted_kernel)

    extinction = receiver_extinction_factor(
        Z, rho, h, snap["unit_length_cgs"], snap["unit_mass_cgs"], sigma_d_band_cgs
    )
    injected_power = float(np.sum(weight * L_band * extinction))
    S_true = weight * L_band * extinction / mass
    return S_true, injected_power, rs


def discrete_steady_state(
    snap, lam, L_band, sigma_d_band_cgs, n_iter, iter_tol, relax, c_hyp_dial=0.1
):
    """Iterate the implemented staggered exact-relaxation update to its own
    fixed point, on the run's REAL particle positions/h/rho/mass, for one
    band's physical screening length `lam`.

    The per-iteration ordering mirrors the C code's own task ordering
    exactly: `div(F)` from the previous iteration's `F` (density loop), then
    the `u` update (radiation_end_density_propagation), then the closure
    tensor built from the NEW `u` and the OLD `F` (gradient loop), then the
    `F` update (radiation_end_gradient_propagation), then the flux limiter
    against the new `u` (radiation_end_force_propagation).

    `c_hyp` and `dt` are set to an arbitrary, numerically convenient value.
    The fixed point is independent of both, under the M1 closure as it was
    under P1: at the fixed point `u = (lam/c_hyp)*(S - div F)` and
    `F = -c_hyp*lam*g`, so the reduced flux `f = |F|/(c_hyp*u) = lam*|g|/u`
    and the limiter condition `|F| <= c_hyp*u` are both `c_hyp`-free, and
    with the source carrying its own `c_hyp/c` rescale the `c_hyp` in
    `u = (lam/c_hyp)*S_used` cancels too, leaving `u = lam*S_true/c`. Run
    the script with `--check-chyp-invariance` to re-measure this rather
    than take it on trust.

    Returns (u, F, f, n_iterations, converged).
    """
    pos, h, rho, mass = snap["pos"], snap["h"], snap["rho"], snap["mass"]
    boxsize = snap["boxsize"]
    N = len(pos)
    ii, jj, dx, r, wi_dr, wj_dr = build_pairs(pos, h, boxsize)

    c_hyp = c_hyp_dial * float(np.median(h))
    a = c_hyp / lam  # = c_hyp*kappa*rho*dt, with dt = 1 internal time unit
    e_, ph = np.exp(-a), phi_relaxation_factor(a)

    S_true, _, _ = injection_source(snap, sigma_d_band_cgs, L_band)
    c_light_internal = C_LIGHT_CGS * snap["unit_time_cgs"] / snap["unit_length_cgs"]
    # The c_hyp/c source rescale, applied exactly where the C code applies
    # it (radiation_end_density_propagation), with this solve's own c_hyp.
    S = S_true * (c_hyp / c_light_internal)

    u = np.zeros(N)
    F = np.zeros((3, N))
    it = 0
    resid = np.inf
    for it in range(n_iter):
        dF = div_F(F, ii, jj, dx, r, wi_dr, wj_dr, rho, mass)
        u_new = e_ * u + ph * (S - dF)
        D, _ = m1_closure_tensor(u_new, F, c_hyp)
        g = grad_P(u_new, D, ii, jj, dx, r, wi_dr, wj_dr, rho, mass)
        F_new = e_ * F - ph * c_hyp**2 * g
        F_new = apply_flux_limiter(u_new, F_new, c_hyp)

        if relax != 1.0:
            u_new = u + relax * (u_new - u)
            F_new = F + relax * (F_new - F)
            F_new = apply_flux_limiter(u_new, F_new, c_hyp)

        # Convergence on BOTH fields: the limiter acts on F, so a converged
        # u does not by itself prove F has converged.
        du = np.max(np.abs(u_new - u)) / max(np.max(np.abs(u_new)), 1e-300)
        dFrel = np.max(np.abs(F_new - F)) / max(np.max(np.abs(F_new)), 1e-300)
        resid = max(du, dFrel)
        u, F = u_new, F_new
        if resid < iter_tol and it > 20:
            break

    _, f = m1_closure_tensor(u, F, c_hyp)
    return u, F, f, it + 1, resid < iter_tol


def main():
    opt = parse_options()
    files = sorted(glob.glob(opt.snapshot))
    if not files:
        raise RuntimeError(f"No snapshots match {opt.snapshot}")

    # Use the last snapshot: the run is designed to reach a converged
    # steady state well before it ends (see README).
    snap = load_snapshot(files[-1])

    r = radial_distance(snap["pos"], snap["star_pos"], snap["boxsize"])
    Z_mean = float(np.mean(snap["Z"]))
    rho_mean = float(np.mean(snap["rho"]))
    h_mean = float(np.median(snap["h"]))

    lambda_fuv = (
        analytic_lambda_cgs(
            Z_mean,
            rho_mean,
            snap["unit_length_cgs"],
            snap["unit_mass_cgs"],
            SIGMA_D_FUV_CGS,
        )
        / snap["unit_length_cgs"]
    )
    lambda_lw = (
        analytic_lambda_cgs(
            Z_mean,
            rho_mean,
            snap["unit_length_cgs"],
            snap["unit_mass_cgs"],
            SIGMA_D_LW_CGS,
        )
        / snap["unit_length_cgs"]
    )

    # Bin edges: exclude the star's own kernel (r too small, where discrete
    # injection geometry dominates) and the outer quarter of the half-box
    # (periodic images start to bias the minimum-image distance there).
    # Identical bins for the simulation and the discrete-solve prediction.
    half_box = snap["boxsize"] / 2.0
    r_min = 2.0 * (snap["boxsize"] / snap["pos"].shape[0] ** (1.0 / 3.0))
    r_max = 0.7 * half_box
    edges = np.linspace(r_min, r_max, opt.n_bins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])

    print(f"Snapshot: {files[-1]} (t={snap['time']:.4e})")
    print(f"Mean Z={Z_mean:.4e}, mean rho={rho_mean:.4e} (internal units)")
    print(f"h/lambda: FUV={h_mean / lambda_fuv:.3f}, LW={h_mean / lambda_lw:.3f}")
    print(f"Fit/compare radial range: [{r_min:.4e}, {r_max:.4e}] (internal units)")
    print()
    print("Solving each band's discrete steady state on the run's own real")
    print("particle positions/h/rho/mass (may take a few seconds)...")

    c_light_internal = C_LIGHT_CGS * snap["unit_time_cgs"] / snap["unit_length_cgs"]
    results = {}
    for band, lam, sigma_d, L_band, u_sim in (
        ("FUV", lambda_fuv, SIGMA_D_FUV_CGS, snap["L_FUV"], snap["u_fuv"]),
        ("LW", lambda_lw, SIGMA_D_LW_CGS, snap["L_LW"], snap["u_lw"]),
    ):
        u_pred, _, f_pred, iters, converged = discrete_steady_state(
            snap, lam, L_band, sigma_d, opt.max_iter, opt.iter_tol, opt.relax
        )
        _, injected_power, _ = injection_source(snap, sigma_d, L_band)
        results[band] = dict(
            lam=lam,
            u_sim=u_sim,
            u_pred=u_pred,
            f_pred=f_pred,
            iters=iters,
            converged=converged,
            injected_power=injected_power,
            binned_sim=radial_bin(r, u_sim, edges),
            binned_pred=radial_bin(r, u_pred, edges),
        )

    ok = True
    for band, res in results.items():
        lam = res["lam"]
        print()
        if not res["converged"]:
            print(
                f"{band}: WARNING -- the discrete fixed-point iteration did not "
                f"converge in {res['iters']} iterations; the prediction below is "
                f"unreliable. Try a smaller --relax (see its help text)."
            )
            ok = False

        # --- Gate 1: the binned profile, simulation vs. discrete prediction.
        valid = (
            ~np.isnan(res["binned_sim"])
            & ~np.isnan(res["binned_pred"])
            & (res["binned_pred"] > 0)
        )
        if valid.sum() < 3:
            print(f"{band}: too few valid radial bins to compare.")
            ok = False
            continue
        rel = (
            np.abs(res["binned_sim"][valid] - res["binned_pred"][valid])
            / res["binned_pred"][valid]
        )
        med_rel, max_rel = float(np.median(rel)), float(np.max(rel))
        profile_ok = med_rel < opt.tol
        ok = ok and profile_ok
        print(
            f"{band}: profile u(r) vs. discrete prediction over {valid.sum()} bins: "
            f"median rel_err={med_rel:.4f}, max={max_rel:.4f} "
            f"-> {'PASS' if profile_ok else 'FAIL'} (tol={opt.tol}) "
            f"[discrete solve: {res['iters']} iterations, "
            f"{'converged' if res['converged'] else 'NOT converged'}]"
        )

        # --- Gate 2: the closure-independent amplitude identity.
        lam_part = (
            analytic_lambda_cgs(
                snap["Z"],
                snap["rho"],
                snap["unit_length_cgs"],
                snap["unit_mass_cgs"],
                SIGMA_D_FUV_CGS if band == "FUV" else SIGMA_D_LW_CGS,
            )
            / snap["unit_length_cgs"]
        )
        lhs_sim = float(np.sum(snap["mass"] * res["u_sim"] / lam_part))
        lhs_pred = float(np.sum(snap["mass"] * res["u_pred"] / lam_part))
        rhs = res["injected_power"] / c_light_internal
        amp_rel_sim = abs(lhs_sim - rhs) / rhs
        amp_rel_pred = abs(lhs_pred - rhs) / rhs
        amp_ok = amp_rel_sim < opt.amp_tol
        ok = ok and amp_ok
        print(
            f"{band}: amplitude identity sum(m*u/lambda) = (1/c)*sum(w*L*ext): "
            f"sim={lhs_sim:.4e}, discrete={lhs_pred:.4e}, analytic={rhs:.4e}, "
            f"rel_err(sim)={amp_rel_sim:.4f}, rel_err(discrete)={amp_rel_pred:.4f} "
            f"-> {'PASS' if amp_ok else 'FAIL'} (tol={opt.amp_tol})"
        )

        # --- Informational: which closure branch the run is actually on,
        # and the old lambda-fit numbers, for continuity with earlier logs.
        in_range = (r > r_min) & (r < r_max)
        f_range = res["f_pred"][in_range]
        print(
            f"{band}: reduced flux f in the compare range (discrete solve): "
            f"median={np.median(f_range):.3f}, "
            f"p10={np.percentile(f_range, 10):.3f}, "
            f"p90={np.percentile(f_range, 90):.3f} "
            f"(f -> 1 is free streaming, f -> 0 is the isotropic/diffusive branch)"
        )
        _, lam_sim_fit = fit_slope(
            centres[valid], res["binned_sim"][valid], r_min, r_max
        )
        _, lam_pred_fit = fit_slope(
            centres[valid], res["binned_pred"][valid], r_min, r_max
        )
        print(
            f"{band}: (informational, NOT gated -- see this script's docstring) "
            f"semi-log fit of u*r: sim={lam_sim_fit:.4e}, "
            f"discrete={lam_pred_fit:.4e}, lambda_analytic={lam:.4e}"
        )

    # Informational only: global negativity and profile monotonicity.
    n_negative = int(np.sum(snap["u_fuv"] < 0) + np.sum(snap["u_lw"] < 0))
    print()
    print(
        f"(informational) particles with u < 0: "
        f"{n_negative}/{2 * len(snap['u_fuv'])}"
    )

    if opt.check_chyp_invariance:
        print()
        print("Re-solving FUV at a 5x smaller arbitrary c_hyp...")
        u_alt, _, _, _, conv_alt = discrete_steady_state(
            snap,
            lambda_fuv,
            snap["L_FUV"],
            SIGMA_D_FUV_CGS,
            opt.max_iter,
            opt.iter_tol,
            opt.relax,
            c_hyp_dial=0.02,
        )
        ref = results["FUV"]["u_pred"]
        scale = np.max(np.abs(ref))
        print(
            f"c_hyp-invariance of the discrete fixed point: "
            f"max|du|/max|u| = {np.max(np.abs(u_alt - ref)) / scale:.3e} "
            f"({'converged' if conv_alt else 'NOT converged'})"
        )

    fig, ax = plt.subplots(figsize=(6, 5))
    for band, colour, marker in (("FUV", "C0", "o"), ("LW", "C1", "s")):
        res = results[band]
        v = ~np.isnan(res["binned_sim"]) & ~np.isnan(res["binned_pred"])
        ax.semilogy(
            centres[v],
            res["binned_sim"][v],
            marker + "-",
            color=colour,
            label=f"{band}: sim",
        )
        ax.semilogy(
            centres[v],
            res["binned_pred"][v],
            "--",
            color=colour,
            label=f"{band}: discrete-solve prediction",
        )
    ax.set_xlabel("r (internal length units)")
    ax.set_ylabel(r"$u(r)$")
    ax.legend(fontsize=7, loc="upper right")
    fig.tight_layout()
    fig.savefig(opt.output, dpi=150)
    print(f"Plot saved to {opt.output}")

    if not ok:
        sys.exit(1)


if __name__ == "__main__":
    main()
