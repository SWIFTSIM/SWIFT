/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published
 * by the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 ******************************************************************************/
#ifndef SWIFT_RADIATION_PROPAGATION_IACT_GEAR_H
#define SWIFT_RADIATION_PROPAGATION_IACT_GEAR_H

/**
 * @file src/feedback/GEAR/radiation_propagation_iact.h
 * @brief Gas-gas density-loop, gradient-loop and force-loop hooks for the
 * hyperbolic M1-relaxation propagation of the per-band u and specific_flux
 * fields.
 *
 * Three pairwise SPH operators are accumulated here, in two different loops,
 * plus the negativity trigger's reference field in a third:
 *
 * - `div(F)` (force loop, accumulated together with the dissipation below in
 *   `runner_iact_[nonsym_]isrf_dissipation`): the shared-coefficient
 *   construction mirroring `src/rt/SPHM1RT/rt_gradients.h`'s
 *   `radiation_divergence_SPH` `diffmode==1` branch. A single shared scalar,
 *   built from both particles' own kernel-gradient terms, is applied with
 *   mirrored mass and sign to each side, so it is exactly mass-conserving
 *   under transport alone for any h_i != h_j, rho_i != rho_j, PROVIDED both
 *   sides of every pair receive it. The density and gradient loops are
 *   type-1 loops: they reach particle i only for r < H_i, so a pair with
 *   H_i <= r < H_j would credit j and never debit i. The force loop is a
 *   type-2 loop and fires both sides whenever either kernel reaches, which
 *   is why the divergence lives there. It reads this step's flux F^{n+1},
 *   already relaxed by the extra ghost. Closure-independent.
 * - `grad(u)` (gradient loop, `runner_iact_[nonsym_]isrf_gradient`): the
 *   anisotropic M1 pressure-tensor divergence, `diffmode==0` form (each
 *   particle's own separate `wi_dr`/`wj_dr`, no shared average, and no
 *   grad-h `forcef` factor: the divergence loop above carries none
 *   either, which the pairing below requires, unlike
 *   `radiation_gradient_aniso_SPH`'s `diffmode==2` branch in
 *   `src/rt/SPHM1RT/rt_gradients.h`/`rt_iact.h:582-627`, which SPHM1RT uses
 *   instead). This own-derivative form is the COMPLEMENT of the
 *   divergence loop's shared-coefficient (`diffmode==1`) construction
 *   above, not a copy of it: `diffmode==1` divergence paired with
 *   `diffmode==0` gradient is what makes the two exactly skew-adjoint
 *   (see below). Matching both loops to the same `diffmode` would NOT
 *   achieve this. Per particle D(f) tensor
 *   (#radiation_get_m1_closure_tensor_band). This is the exact skew-adjoint
 *   of the `diffmode==1` divergence above, in the `D^-1`-weighted inner
 *   product (not the plain `m*rho` one the P1 pairing used), for any
 *   `h_i != h_j`, `rho_i != rho_j`, whenever `D` is locally constant
 *   between neighbours (`D_i = D_j`), which includes `F=0`, every
 *   particle's initial condition and permanent far-field state; a residual
 *   `O(h*grad(D))` remains when `D` varies between neighbours. This removes
 *   only that kernel-derivative-mismatch residual `R_h`. In the plain
 *   isotropic metric, a separate, irremovable anisotropy residual `R_D`
 *   remains regardless of `h` or of whether `D_i = D_j`: nonzero whenever
 *   the closure is anisotropic at all, including at uniform `h` with the
 *   same constant anisotropic `D` on both particles. `R_D` is a PDE-level
 *   property of M1 itself (the isotropic-metric functional this fix
 *   generalizes is not a conserved quantity of M1 for anisotropic `D`),
 *   not a discretization defect, and this fix does not address it.
 *   Whether the staggered exact-relaxation time integrator's stability
 *   argument needs more than this weaker adjointness remains an open
 *   question, not resolved here. Unlike the divergence, this operator is
 *   complete per particle under a type-1 loop: i's own term carries `wi_dr`,
 *   which vanishes for r >= H_i, and nothing is mirrored onto j.
 * - The negativity-triggered artificial dissipation (force loop,
 *   `runner_iact_[nonsym_]isrf_dissipation`): a triggered pairwise
 *   conductivity on the `rho*u` jump, credited to one particle and debited
 *   from the other. It lives in the force loop for the same pair-coverage
 *   reason as the divergence, and shares that hook's kernel evaluations.
 *   SPHENIX's own artificial viscosity lives in the force loop for the
 *   identical reason. Closure-independent.
 * - The trigger's reference field, the kernel mean of the neighbours'
 *   `|rho*u|` (density loop, `runner_iact_[nonsym_]isrf_propagation`): a
 *   per-particle mean with no mirrored term, so the type-1 dispatch is
 *   complete for it.
 *
 * Time levels within one step. Nothing writes `u` between the drift snapshot
 * and the end-force ghost (injection with propagation on writes only the
 * dose reservoir), so the live `u` the gradient loop and the dissipation
 * read IS `u^n`, equal to `u_prev`. The order is: the gradient loop builds
 * `grad(u^n)` with `D(u^n, F^n)`; the extra ghost relaxes `F^{n+1}` and
 * limits it against `u^n`; the force loop accumulates `div(F^{n+1})` and the
 * dissipation of `u^n`; the end-force ghost produces `u^{n+1}`.
 *
 * Comoving-to-physical convention. Every quantity SWIFT hands these loops is
 * comoving (`dx`, `r`, `h`, and the `rho_prev` snapshot), while every field
 * they accumulate into is PHYSICAL, like `u`, the tracked flux, `c_hyp` and
 * `kappa`. The three spatial operators above each carry exactly one net
 * inverse length: with comoving length `= physical/a`, comoving density
 * `= physical*a^dim` and a kernel derivative scaling as `h^-(dim+1)`,
 * evaluating any of them on comoving inputs returns `a` times the physical
 * value, in any dimension (divergence: `dx` gives `a^-1`, `1/rho` gives
 * `a^-dim`, `wi_dr` gives `a^(dim+1)` and `r_inv` gives `a^1`, summing to
 * `a^1`; the gradient and the negativity-triggered dissipation give the same
 * total term
 * by term). Each is therefore closed with a single named conversion factor,
 * `a_factor_comoving_to_physical = 1/a`, computed once per pair dispatch in
 * the hooks below and applied where the comoving estimate becomes the
 * physical accumulator value, in the manner of `cosmology.c`'s `a_factor_*`
 * scalars. One operator here needs no conversion at all, for the reason
 * given at its own definition: #radiation_dissipation_reference_accumulate_band
 * (read only as a ratio).
 *
 * Every operator reads `p->feedback_data.rho_prev`, a comoving density
 * snapshot cached once per step by `radiation_snapshot_part_propagation`,
 * the same snapshot in every loop: the density loop runs interleaved with
 * SPH's own density accumulation, where `p->rho` is a partial sum, and the
 * skew-adjoint pairing of the divergence and the gradient needs both built
 * from the same `rho_i`/`rho_j`.
 */

#include "dimension.h"
#include "kernel_hydro.h"
#include "radiation.h"

#include <math.h>

/**
 * @brief Band-specific pairwise contribution to particle i's `div(F)`
 * accumulator, and mirrored (mass-weighted, opposite sign) contribution to
 * particle j's, from a single shared coefficient built from both particles'
 * own flux, density, and kernel-gradient terms.
 *
 * @param dx Comoving separation vector (pi - pj).
 * @param r_inv Inverse comoving particle separation.
 * @param wi_dr Particle i's own kernel-gradient term,
 * h_i^-(dim+1) * dW/dq|_{r/h_i}.
 * @param wj_dr Particle j's own kernel-gradient term,
 * h_j^-(dim+1) * dW/dq|_{r/h_j}.
 * @param mi Particle i's mass.
 * @param mj Particle j's mass.
 * @param rho_i Particle i's cached comoving density snapshot.
 * @param rho_j Particle j's cached comoving density snapshot.
 * @param F_i Particle i's tracked flux (this band).
 * @param F_j Particle j's tracked flux (this band).
 * @param a_factor_comoving_to_physical `1/a`, the file header's single
 * conversion factor: the comoving inputs make the shared coefficient `a`
 * times the physical divergence, and both accumulators are physical.
 * @param div_F_i (return, accumulated) Particle i's div(F) accumulator.
 * @param div_F_j (return, accumulated) Particle j's div(F) accumulator.
 */
__attribute__((always_inline)) INLINE static void
radiation_divergence_accumulate_band(const float dx[3], float r_inv,
                                     float wi_dr, float wj_dr, float mi,
                                     float mj, float rho_i, float rho_j,
                                     const float F_i[3], const float F_j[3],
                                     float a_factor_comoving_to_physical,
                                     float *div_F_i, float *div_F_j) {

  const float Fi_dot_dx = F_i[0] * dx[0] + F_i[1] * dx[1] + F_i[2] * dx[2];
  const float Fj_dot_dx = F_j[0] * dx[0] + F_j[1] * dx[1] + F_j[2] * dx[2];

  const float Phi_ij =
      (Fi_dot_dx / rho_i * wi_dr * r_inv + Fj_dot_dx / rho_j * wj_dr * r_inv) *
      a_factor_comoving_to_physical;

  *div_F_i += mj * Phi_ij;
  *div_F_j += -mi * Phi_ij;
}

/**
 * @brief Band-specific pairwise contribution to each particle's kernel-mean
 * `|rho_prev*u_prev|` reference accumulator, the local field scale the
 * negativity trigger divides an undershoot by (radiation_isrf.c's
 * #radiation_update_dissipation_alpha_band).
 *
 * Accumulated in the density loop, from the stable `u_*_prev` snapshot and
 * `rho_prev`, so the value the trigger reads does not drift across a
 * particle's h-iterations. The dissipation term the trigger drives is
 * accumulated separately, in the force loop
 * (#radiation_dissipation_force_accumulate_band).
 *
 * Takes no comoving-to-physical conversion, unlike the three spatial
 * operators in this file: its only consumer divides it into `rho_prev*u`,
 * which carries the same `a^3` comoving-density weight, so the two cancel.
 * Converting here would introduce a bias rather than remove one.
 *
 * @param wi Particle i's own kernel value, W(r/h_i)*h_i^-dim.
 * @param wj Particle j's own kernel value, W(r/h_j)*h_j^-dim.
 * @param mi Particle i's mass.
 * @param mj Particle j's mass.
 * @param rho_i Particle i's cached comoving density snapshot.
 * @param rho_j Particle j's cached comoving density snapshot.
 * @param u_i_prev Particle i's snapshotted specific field (this band).
 * @param u_j_prev Particle j's snapshotted specific field (this band).
 * @param ngb_mean_abs_u_V_i (return, accumulated) Particle i's kernel-mean
 * `|rho_prev*u_prev|` accumulator.
 * @param ngb_mean_abs_u_V_j (return, accumulated) Particle j's kernel-mean
 * `|rho_prev*u_prev|` accumulator.
 */
__attribute__((always_inline)) INLINE static void
radiation_dissipation_reference_accumulate_band(float wi, float wj, float mi,
                                                float mj, float rho_i,
                                                float rho_j, float u_i_prev,
                                                float u_j_prev,
                                                float *ngb_mean_abs_u_V_i,
                                                float *ngb_mean_abs_u_V_j) {

  *ngb_mean_abs_u_V_i += (mj / rho_j) * wi * fabsf(rho_j * u_j_prev);
  *ngb_mean_abs_u_V_j += (mi / rho_i) * wj * fabsf(rho_i * u_i_prev);
}

/**
 * @brief Band-specific pairwise contribution to particle i's
 * negativity-triggered artificial-dissipation source term, and mirrored
 * (mass-weighted, opposite sign) contribution to particle j's.
 *
 * `v_sig,ij = alpha_ij * min(c_hyp_i, c_hyp_j)` is a signal VELOCITY, with
 * no `h` factor: the length scale enters only through `Wbar_ij`'s own
 * `h^-(dim+1)` normalisation, exactly as for an ordinary SPH Laplacian.
 *
 * The pair coefficient combines the two components the extra ghost stores
 * separately (#feedback_isrf_band_data.dissipation_alpha_trigger,
 * #feedback_isrf_band_data.dissipation_alpha_floor):
 *
 *   `alpha_ij = max(trigger_i, trigger_j, floor_i, floor_j)`.
 *
 * Both components enter unconditionally. The trigger is local and reactive:
 * it only ever fires on a particle whose own field has already gone wrong.
 * The floor is anticipatory and carries no per-pair condition of its own;
 * it is gated only per particle, on the resolution ratio `h/lambda`
 * (#radiation_dissipation_alpha_floor_band), so that it supplies
 * dissipation on a positive front the trigger is structurally blind to and
 * switches itself off where the absorption length is resolved.
 *
 * Positivity of an explicit hyperbolic update is a property of the operator
 * and the CFL number, not of the field's local shape, so the dissipation
 * coefficient itself is the wrong place to discriminate between a
 * sub-resolution oscillation and a resolved contrast: a first-order upwind
 * flux with `~v_sig*|d_ij|/2` dissipation is unconditionally positive, and
 * its cost is over-diffusion at sharp, resolved contrasts. That cost is a
 * disclosed limitation of the scheme, shared with the wider M1/P1
 * literature.
 *
 * Because `alpha_ij <= max(alpha_max, alpha_floor)`, the joint
 * `(alpha, C_hyp)` stability bound checked in feedback_props_init() holds.
 *
 * Accumulated in the force loop, after the extra ghost has set this step's
 * `alpha`, from the live `u`, which is still `u^n` there (see the file
 * header). The force loop runs exactly once per step, so there is no
 * h-iteration stability requirement here.
 *
 * The credit/debit pair is applied unconditionally, with no mutual-reach
 * gate: the force loop's dispatch fires both sides whenever either kernel
 * reaches, so the pair is always mirrored and `sum_i m_i*dissipation_u_i`
 * is exactly zero over a same-bin active-active pair for any h_i, h_j.
 * `Psi_ji = -Psi_ij` under the exchange `i <-> j`, `dx -> -dx`: `Wbar_ij`
 * is symmetric and `d_ij` is antisymmetric.
 *
 * @param wi_dr See #radiation_divergence_accumulate_band.
 * @param wj_dr See #radiation_divergence_accumulate_band.
 * @param mi Particle i's mass.
 * @param mj Particle j's mass.
 * @param rho_i Particle i's cached comoving density snapshot.
 * @param rho_j Particle j's cached comoving density snapshot.
 * @param c_hyp_i Particle i's own hyperbolic propagation speed.
 * @param c_hyp_j Particle j's own hyperbolic propagation speed.
 * @param alpha_trigger_i Particle i's
 * #feedback_isrf_band_data.dissipation_alpha_trigger (this band).
 * @param alpha_trigger_j Particle j's
 * #feedback_isrf_band_data.dissipation_alpha_trigger (this band).
 * @param alpha_floor_i Particle i's
 * #feedback_isrf_band_data.dissipation_alpha_floor (this band).
 * @param alpha_floor_j Particle j's
 * #feedback_isrf_band_data.dissipation_alpha_floor (this band).
 * @param u_i Particle i's live specific field `u^n` (this band).
 * @param u_j Particle j's live specific field `u^n` (this band).
 * @param a_factor_comoving_to_physical `1/a`, the file header's single
 * conversion factor, applied to the shared coefficient `Psi_ij`.
 * @param dissipation_u_i (return, accumulated) Particle i's dissipation
 * source-term accumulator.
 * @param dissipation_u_j (return, accumulated) Particle j's dissipation
 * source-term accumulator.
 */
__attribute__((always_inline)) INLINE static void
radiation_dissipation_force_accumulate_band(
    float wi_dr, float wj_dr, float mi, float mj, float rho_i, float rho_j,
    float c_hyp_i, float c_hyp_j, float alpha_trigger_i, float alpha_trigger_j,
    float alpha_floor_i, float alpha_floor_j, float u_i, float u_j,
    float a_factor_comoving_to_physical, float *dissipation_u_i,
    float *dissipation_u_j) {

  const float d_ij = rho_i * u_i - rho_j * u_j;

  /* Split out rather than nested: max() expands to a statement expression
   * with its own locals, which -Wshadow rejects when nested. */
  const float alpha_trigger_ij = max(alpha_trigger_i, alpha_trigger_j);
  const float alpha_floor_ij = max(alpha_floor_i, alpha_floor_j);
  const float alpha_ij = max(alpha_trigger_ij, alpha_floor_ij);

  const float Wbar_ij = 0.5f * (wi_dr + wj_dr);
  const float v_sig_ij = alpha_ij * min(c_hyp_i, c_hyp_j);
  const float Psi_ij = v_sig_ij * d_ij * Wbar_ij / (rho_i * rho_j) *
                       a_factor_comoving_to_physical;

  *dissipation_u_i += mj * Psi_ij;
  *dissipation_u_j += -mi * Psi_ij;
}

/**
 * @brief M1 closure tensor `D(f)` for one particle, one band, built from its
 * own `(u, F, c_M)`. `c_M` is the
 * same speed already carried as #feedback_part_data.c_hyp, reinterpreted
 * as the fastest M1 characteristic, `f=1`, not a new field.
 *
 * `f = min(1, |F|/(c_M*u))` for `u > 0`, `f = 0` for `u <= 0`;
 * `chi(f) = (3+4f^2)/(5+2*sqrt(4-3f^2))`;
 * `D(f) = (1-chi)/2 I + (3chi-1)/2 (n dyadic n)`, `n = F/|F|`.
 *
 * Zero-flux guard, mandatory: `F = 0` is every particle's initial condition
 * and permanent far-field state, not a corner case. `F2 = F.F`,
 * `F_inv = (F2 > 0) ? 1/sqrt(F2) : 0`, `n = F*F_inv`: a zero-guarded
 * reciprocal, so `n` stays well-defined at `F = 0` instead of dividing by
 * zero. `f`'s own division is guarded the same way: `c_M*u` is computed once
 * and only divided into when it is strictly positive, which also folds in the
 * `u <= 0` case (`f = 0`) without a separate branch. At `F = 0`, `f = 0`,
 * `chi = 1/3`, the `(3*chi-1)/2 = 0` coefficient multiplies the guarded,
 * well-defined zero `n` rather than a NaN.
 *
 * `F2`, `|F|`, `c_M*u` and `f` are formed in double: in float32, `F.F`
 * underflows to zero once `|F| < sqrt(FLT_MIN) ~ 1.1e-19` (internal units),
 * which would take the zero-flux branch for a nonzero flux and turn a beam
 * into an isotropic closure. Any nonzero float32 `F` has a nonzero double
 * `F2`, so `F2 > 0` in double is exactly `F != 0`. `n` and the closure
 * coefficients are O(1) and stay in float.
 *
 * @param u This band's specific field `u^n` for this particle.
 * @param F This particle's tracked flux (this band).
 * @param c_M This particle's own #feedback_part_data.c_hyp.
 * @param D (return) The 3x3 closure tensor.
 */
__attribute__((always_inline)) INLINE static void
radiation_get_m1_closure_tensor_band(float u, const float F[3], float c_M,
                                     float D[3][3]) {

  const double F2 = (double)F[0] * (double)F[0] + (double)F[1] * (double)F[1] +
                    (double)F[2] * (double)F[2];
  const double F_inv = (F2 > 0.) ? 1. / sqrt(F2) : 0.;
  const double Fmag = F2 * F_inv; /* sqrt(F2), no second sqrt call */
  const float n[3] = {(float)(F[0] * F_inv), (float)(F[1] * F_inv),
                      (float)(F[2] * F_inv)};

  const double denom = (double)c_M * (double)u;
  const float f = (denom > 0.) ? (float)min(Fmag / denom, 1.) : 0.f;

  const float sq = 4.f - 3.f * f * f;
  const float chi = (3.f + 4.f * f * f) / (5.f + 2.f * sqrtf(sq));

  const float iso_coeff = 0.5f * (1.f - chi);
  const float aniso_coeff = 0.5f * (3.f * chi - 1.f);

  for (int a = 0; a < 3; a++) {
    D[a][0] = aniso_coeff * n[a] * n[0];
    D[a][1] = aniso_coeff * n[a] * n[1];
    D[a][2] = aniso_coeff * n[a] * n[2];
    D[a][a] += iso_coeff;
  }
}

/**
 * @brief Band-specific pairwise contribution to particle i's `grad(u)`
 * accumulator (and mirrored contribution to particle j's), the anisotropic
 * M1 pressure-tensor divergence `1/rho * div(D(f)*rho*u)`.
 *
 * `diffmode == 0` form (each particle's own separate `wi_dr`/`wj_dr`, no
 * shared average, no grad-h `forcef` factor): `tempi - tempj` on
 * `D_i`/`D_j` applied with each particle's own kernel derivative, matching
 * this project's existing `d_ij = rho_i*u_i - rho_j*u_j` jump construction
 * in spirit (own-derivative, no averaging) but NOT matching
 * `radiation_divergence_accumulate_band` above, which uses `diffmode==1`
 * (a single shared coefficient, not each particle's own derivative) --
 * the two loops are deliberately COMPLEMENTARY, not identical: pairing
 * `diffmode==1` divergence with `diffmode==0` gradient is what makes them
 * exactly skew-adjoint (see below); using the same `diffmode` for both
 * would not. Also not `src/rt/SPHM1RT/rt_gradients.h`'s
 * `radiation_gradient_aniso_SPH` `diffmode==2` branch
 * (`src/rt/SPHM1RT/rt_iact.h:582-627`), which uses a shared averaged
 * derivative and (in its own `diffmode==0` branch) a grad-h correction
 * this operator deliberately omits: see the file header comment for why
 * (exact skew-adjointness with the divergence loop, in the `D^-1` metric,
 * whenever `D` is locally constant). `D_i`, `D_j` reduce to
 * `(1/3) I` at `f=0` (both particles' fluxes zero, the isotropic P1
 * limit), so this reduces to the old scalar form's structure with the
 * `1/3` now explicit rather than folded away.
 *
 * @param dx Comoving separation vector (pi - pj).
 * @param r_inv Inverse comoving particle separation.
 * @param wi_dr See #radiation_divergence_accumulate_band.
 * @param wj_dr See #radiation_divergence_accumulate_band.
 * @param mi Particle i's mass.
 * @param mj Particle j's mass.
 * @param rho_i Particle i's cached comoving density snapshot.
 * @param rho_j Particle j's cached comoving density snapshot.
 * @param u_i Particle i's specific field `u^n` (this band).
 * @param u_j Particle j's specific field `u^n` (this band).
 * @param D_i Particle i's own M1 closure tensor (this band), from
 * #radiation_get_m1_closure_tensor_band.
 * @param D_j Particle j's own M1 closure tensor (this band).
 * @param a_factor_comoving_to_physical `1/a`, the file header's single
 * conversion factor: folded into `fac_i`/`fac_j` so both accumulators come
 * out physical.
 * @param grad_u_i (return, accumulated) Particle i's grad(u) accumulator.
 * @param grad_u_j (return, accumulated) Particle j's grad(u) accumulator.
 */
__attribute__((always_inline)) INLINE static void
radiation_gradient_accumulate_band(const float dx[3], float r_inv, float wi_dr,
                                   float wj_dr, float mi, float mj, float rho_i,
                                   float rho_j, float u_i, float u_j,
                                   const float D_i[3][3], const float D_j[3][3],
                                   float a_factor_comoving_to_physical,
                                   float grad_u_i[3], float grad_u_j[3]) {

  const float rho_i_inv = 1.f / rho_i;
  const float rho_j_inv = 1.f / rho_j;

  float temp_i[3], temp_j[3];
  for (int k = 0; k < 3; k++) {
    const float Di_dot_dx =
        D_i[k][0] * dx[0] + D_i[k][1] * dx[1] + D_i[k][2] * dx[2];
    const float Dj_dot_dx =
        D_j[k][0] * dx[0] + D_j[k][1] * dx[1] + D_j[k][2] * dx[2];
    temp_i[k] = Di_dot_dx * rho_i * u_i * r_inv;
    temp_j[k] = Dj_dot_dx * rho_j * u_j * r_inv;
  }

  /* Own kernel derivative per particle, no shared average and no grad-h
   * `forcef` factor: restores exact adjointness with the divergence loop
   * above (D^-1 metric, D locally constant). See the header comment. */
  const float fac_i =
      mj * rho_i_inv * rho_i_inv * wi_dr * a_factor_comoving_to_physical;
  const float fac_j =
      mi * rho_j_inv * rho_j_inv * wj_dr * a_factor_comoving_to_physical;

  for (int k = 0; k < 3; k++) {
    grad_u_i[k] += -(temp_i[k] - temp_j[k]) * fac_i;
    grad_u_j[k] += -(temp_i[k] - temp_j[k]) * fac_j;
  }
}

/**
 * @brief Density-loop propagation interaction between two particles
 * (symmetric): both particles' trigger reference accumulators are updated.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor (unused, for the same reason as `us` below).
 * @param H Current Hubble parameter (unused, for the same reason as `us`
 * below).
 * @param us Unit system (unused: the kernel mean needs only positions,
 * masses, the cached density snapshot, and the snapshotted field).
 */
__attribute__((always_inline)) INLINE static void runner_iact_isrf_propagation(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const struct unit_system *us) {

  const float r = sqrtf(r2);
  const float hi_inv = 1.f / hi;
  const float hj_inv = 1.f / hj;

  float wi, wj;
  kernel_eval(r * hi_inv, &wi);
  kernel_eval(r * hj_inv, &wj);
  wi *= pow_dimension(hi_inv);
  wj *= pow_dimension(hj_inv);

  struct feedback_part_data *fdi = &pi->feedback_data;
  struct feedback_part_data *fdj = &pj->feedback_data;
  const float rho_i = fdi->rho_prev;
  const float rho_j = fdj->rho_prev;
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  for (int b = 0; b < ISRF_BAND_COUNT; b++) {
    struct feedback_isrf_band_data *bi = &fdi->isrf_band[b];
    struct feedback_isrf_band_data *bj = &fdj->isrf_band[b];

    radiation_dissipation_reference_accumulate_band(
        wi, wj, mi, mj, rho_i, rho_j, bi->u_prev, bj->u_prev,
        &bi->ngb_mean_abs_u_V, &bj->ngb_mean_abs_u_V);
  }
}

/**
 * @brief Density-loop propagation interaction between two particles
 * (non-symmetric): only particle i's trigger reference accumulator is
 * updated.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (its own accumulators not updated).
 * @param a Current scale factor (unused, see #runner_iact_isrf_propagation).
 * @param H Current Hubble parameter (unused, see
 * #runner_iact_isrf_propagation).
 * @param us Unit system (unused, see #runner_iact_isrf_propagation).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_isrf_propagation(const float r2, const float dx[3],
                                    const float hi, const float hj,
                                    struct part *restrict pi,
                                    const struct part *restrict pj,
                                    const float a, const float H,
                                    const struct unit_system *us) {

  const float r = sqrtf(r2);
  const float hi_inv = 1.f / hi;
  const float hj_inv = 1.f / hj;

  float wi, wj;
  kernel_eval(r * hi_inv, &wi);
  kernel_eval(r * hj_inv, &wj);
  wi *= pow_dimension(hi_inv);
  wj *= pow_dimension(hj_inv);

  struct feedback_part_data *fdi = &pi->feedback_data;
  const struct feedback_part_data *fdj = &pj->feedback_data;
  const float rho_i = fdi->rho_prev;
  const float rho_j = fdj->rho_prev;
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  for (int b = 0; b < ISRF_BAND_COUNT; b++) {
    struct feedback_isrf_band_data *bi = &fdi->isrf_band[b];
    const struct feedback_isrf_band_data *bj = &fdj->isrf_band[b];

    /* Particle j's own accumulator is not touched (non-symmetric): pass a
     * discarded local, seeded to 0 rather than read from fdj, as the
     * required (return, accumulated) output. */
    float unused_ngb_mean_abs_u_V = 0.f;

    radiation_dissipation_reference_accumulate_band(
        wi, wj, mi, mj, rho_i, rho_j, bi->u_prev, bj->u_prev,
        &bi->ngb_mean_abs_u_V, &unused_ngb_mean_abs_u_V);
  }
}

/**
 * @brief `grad(u)` interaction between two particles (symmetric): both
 * particles' accumulators are updated.
 *
 * Runs in the gradient loop and reads the live `u`, which is `u^n` there
 * (see the file header): this step's `u` update happens later, in the
 * end-force ghost (radiation_isrf.c's #radiation_end_force_propagation),
 * from the flux the extra ghost relaxes out of this gradient.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_isrf_gradient(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.f / r : 0.f;
  const float hi_inv = 1.f / hi;
  const float hj_inv = 1.f / hj;

  float wi, wi_dx, wj, wj_dx;
  kernel_deval(r * hi_inv, &wi, &wi_dx);
  kernel_deval(r * hj_inv, &wj, &wj_dx);
  const float wi_dr = wi_dx * pow_dimension_plus_one(hi_inv);
  const float wj_dr = wj_dx * pow_dimension_plus_one(hj_inv);
  (void)wi;
  (void)wj;

  struct feedback_part_data *fdi = &pi->feedback_data;
  struct feedback_part_data *fdj = &pj->feedback_data;
  const float rho_i = fdi->rho_prev;
  const float rho_j = fdj->rho_prev;
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  /* Single named conversion factor for every spatial operator below, built
   * once per pair dispatch: see this file's header. */
  const float a_factor_comoving_to_physical = 1.f / a;

  for (int b = 0; b < ISRF_BAND_COUNT; b++) {
    struct feedback_isrf_band_data *bi = &fdi->isrf_band[b];
    struct feedback_isrf_band_data *bj = &fdj->isrf_band[b];

    float D_i[3][3], D_j[3][3];
    radiation_get_m1_closure_tensor_band(bi->u, bi->specific_flux, fdi->c_hyp,
                                         D_i);
    radiation_get_m1_closure_tensor_band(bj->u, bj->specific_flux, fdj->c_hyp,
                                         D_j);

    radiation_gradient_accumulate_band(
        dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, bi->u, bj->u, D_i, D_j,
        a_factor_comoving_to_physical, bi->grad_u, bj->grad_u);
  }
}

/**
 * @brief `grad(u)` interaction between two particles (non-symmetric):
 * only particle i's accumulators are updated.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (its own accumulators not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_isrf_gradient(const float r2, const float dx[3],
                                 const float hi, const float hj,
                                 struct part *restrict pi,
                                 struct part *restrict pj, const float a,
                                 const float H) {

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.f / r : 0.f;
  const float hi_inv = 1.f / hi;
  const float hj_inv = 1.f / hj;

  float wi, wi_dx, wj, wj_dx;
  kernel_deval(r * hi_inv, &wi, &wi_dx);
  kernel_deval(r * hj_inv, &wj, &wj_dx);
  const float wi_dr = wi_dx * pow_dimension_plus_one(hi_inv);
  const float wj_dr = wj_dx * pow_dimension_plus_one(hj_inv);
  (void)wi;
  (void)wj;

  struct feedback_part_data *fdi = &pi->feedback_data;
  const struct feedback_part_data *fdj = &pj->feedback_data;
  const float rho_i = fdi->rho_prev;
  const float rho_j = fdj->rho_prev;
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  /* Single named conversion factor for every spatial operator below, built
   * once per pair dispatch: see this file's header. */
  const float a_factor_comoving_to_physical = 1.f / a;

  for (int b = 0; b < ISRF_BAND_COUNT; b++) {
    struct feedback_isrf_band_data *bi = &fdi->isrf_band[b];
    const struct feedback_isrf_band_data *bj = &fdj->isrf_band[b];

    /* Particle j is `const` here (non-symmetric): its own accumulator is not
     * touched, so pass a discarded, zero-seeded local as the writable
     * destination the shared accumulator function requires for j. */
    float unused_grad_u[3] = {0.f, 0.f, 0.f};

    float D_i[3][3], D_j[3][3];
    radiation_get_m1_closure_tensor_band(bi->u, bi->specific_flux, fdi->c_hyp,
                                         D_i);
    radiation_get_m1_closure_tensor_band(bj->u, bj->specific_flux, fdj->c_hyp,
                                         D_j);

    radiation_gradient_accumulate_band(
        dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, bi->u, bj->u, D_i, D_j,
        a_factor_comoving_to_physical, bi->grad_u, unused_grad_u);
  }
}

/**
 * @brief Force-loop propagation interaction between two particles
 * (symmetric): both particles' `div(F)` and dissipation accumulators are
 * updated.
 *
 * Runs in the force loop, after the extra ghost has relaxed this step's
 * #feedback_isrf_band_data.specific_flux and set
 * #feedback_isrf_band_data.dissipation_alpha_trigger and
 * #feedback_isrf_band_data.dissipation_alpha_floor. The divergence reads that
 * new flux; the dissipation reads the live `u`, still `u^n` (see the file
 * header). The force loop's dispatch fires both sides of a pair whenever
 * either kernel reaches, which is what keeps both mirrored pairs whole at
 * h_i != h_j.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void runner_iact_isrf_dissipation(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H) {

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.f / r : 0.f;
  const float hi_inv = 1.f / hi;
  const float hj_inv = 1.f / hj;

  float wi, wi_dx, wj, wj_dx;
  kernel_deval(r * hi_inv, &wi, &wi_dx);
  kernel_deval(r * hj_inv, &wj, &wj_dx);
  const float wi_dr = wi_dx * pow_dimension_plus_one(hi_inv);
  const float wj_dr = wj_dx * pow_dimension_plus_one(hj_inv);
  (void)wi;
  (void)wj;

  struct feedback_part_data *fdi = &pi->feedback_data;
  struct feedback_part_data *fdj = &pj->feedback_data;
  const float rho_i = fdi->rho_prev;
  const float rho_j = fdj->rho_prev;
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  /* Single named conversion factor for every spatial operator below, built
   * once per pair dispatch: see this file's header. */
  const float a_factor_comoving_to_physical = 1.f / a;

  for (int b = 0; b < ISRF_BAND_COUNT; b++) {
    struct feedback_isrf_band_data *bi = &fdi->isrf_band[b];
    struct feedback_isrf_band_data *bj = &fdj->isrf_band[b];

    radiation_divergence_accumulate_band(
        dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, bi->specific_flux,
        bj->specific_flux, a_factor_comoving_to_physical,
        &bi->div_specific_flux, &bj->div_specific_flux);

    radiation_dissipation_force_accumulate_band(
        wi_dr, wj_dr, mi, mj, rho_i, rho_j, fdi->c_hyp, fdj->c_hyp,
        bi->dissipation_alpha_trigger, bj->dissipation_alpha_trigger,
        bi->dissipation_alpha_floor, bj->dissipation_alpha_floor, bi->u, bj->u,
        a_factor_comoving_to_physical, &bi->dissipation_u, &bj->dissipation_u);
  }
}

/**
 * @brief Force-loop propagation interaction between two particles
 * (non-symmetric): only particle i's `div(F)` and dissipation accumulators
 * are updated.
 *
 * The force loop reaches this variant once per side, so a pair whose two
 * sides are both dispatched still receives both mirrored pairs in full; see
 * #runner_iact_isrf_dissipation.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (its own accumulator not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_isrf_dissipation(const float r2, const float dx[3],
                                    const float hi, const float hj,
                                    struct part *restrict pi,
                                    const struct part *restrict pj,
                                    const float a, const float H) {

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.f / r : 0.f;
  const float hi_inv = 1.f / hi;
  const float hj_inv = 1.f / hj;

  float wi, wi_dx, wj, wj_dx;
  kernel_deval(r * hi_inv, &wi, &wi_dx);
  kernel_deval(r * hj_inv, &wj, &wj_dx);
  const float wi_dr = wi_dx * pow_dimension_plus_one(hi_inv);
  const float wj_dr = wj_dx * pow_dimension_plus_one(hj_inv);
  (void)wi;
  (void)wj;

  struct feedback_part_data *fdi = &pi->feedback_data;
  const struct feedback_part_data *fdj = &pj->feedback_data;
  const float rho_i = fdi->rho_prev;
  const float rho_j = fdj->rho_prev;
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  /* Single named conversion factor for every spatial operator below, built
   * once per pair dispatch: see this file's header. */
  const float a_factor_comoving_to_physical = 1.f / a;

  for (int b = 0; b < ISRF_BAND_COUNT; b++) {
    struct feedback_isrf_band_data *bi = &fdi->isrf_band[b];
    const struct feedback_isrf_band_data *bj = &fdj->isrf_band[b];

    /* Particle j's own accumulators are not touched (non-symmetric): pass
     * discarded locals, seeded to 0 rather than read from fdj, as the
     * required (return, accumulated) outputs. */
    float unused_div_specific_flux = 0.f;
    float unused_dissipation_u = 0.f;

    radiation_divergence_accumulate_band(
        dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, bi->specific_flux,
        bj->specific_flux, a_factor_comoving_to_physical,
        &bi->div_specific_flux, &unused_div_specific_flux);

    radiation_dissipation_force_accumulate_band(
        wi_dr, wj_dr, mi, mj, rho_i, rho_j, fdi->c_hyp, fdj->c_hyp,
        bi->dissipation_alpha_trigger, bj->dissipation_alpha_trigger,
        bi->dissipation_alpha_floor, bj->dissipation_alpha_floor, bi->u, bj->u,
        a_factor_comoving_to_physical, &bi->dissipation_u,
        &unused_dissipation_u);
  }
}

#endif /* SWIFT_RADIATION_PROPAGATION_IACT_GEAR_H */
