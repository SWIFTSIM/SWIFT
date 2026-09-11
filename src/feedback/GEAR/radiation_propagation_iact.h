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
 * hyperbolic M1-relaxation propagation of the u_FUV/specific_flux_FUV (and
 * u_LW/specific_flux_LW) fields.
 *
 * Three pairwise SPH operators are accumulated here, in three different
 * loops:
 *
 * - `div(F)` (density loop, this file's
 * `runner_iact_[nonsym_]isrf_propagation`): the shared-coefficient construction
 * mirroring `src/rt/SPHM1RT/rt_gradients.h`'s `radiation_divergence_SPH`
 *   `diffmode==1` branch. Exactly mass-conserving under transport alone for
 *   any h_i != h_j, rho_i != rho_j (a single shared scalar built from both
 *   particles' own kernel-gradient terms, applied with mirrored mass/sign to
 *   each side). Closure-independent: unchanged by the P1-to-M1 upgrade.
 * - `grad(u)` (gradient loop, `runner_iact_[nonsym_]isrf_gradient`): the
 *   anisotropic M1 pressure-tensor divergence, `diffmode==2` form
 *   (`radiation_gradient_aniso_SPH`'s `diffmode==2` branch,
 *   `src/rt/SPHM1RT/rt_gradients.h`/`rt_iact.h:582-627`): `tempi - tempj` on
 *   a shared, averaged kernel derivative `(wi_dr + wj_dr)*0.5`, per particle
 *   D(f) tensor (#radiation_get_m1_closure_tensor_band). Under P1
 *   (design-lw-fuv-m1-upgrade.md's predecessor), this operator instead used
 *   each particle's own separate `wi_dr`/`wj_dr` (the `diffmode==0` shape)
 *   and was the exact skew-adjoint of the `diffmode==1` divergence above in
 *   the m*rho inner product; the M1 anisotropic form gives up that exact
 *   adjointness (the divergence loop still uses per-particle kernel terms,
 *   this loop no longer does) in exchange for the correct M1 pressure
 *   tensor. Whether the staggered exact-relaxation time integrator's
 *   stability argument still needs that adjointness, or tolerates its loss,
 *   is an open question for Phase 1's stability re-verification, not
 *   resolved here.
 * - The Stage-1 artificial dissipation (force loop,
 *   `runner_iact_[nonsym_]isrf_dissipation`): a triggered pairwise
 *   conductivity on the `rho*u` jump, credited to one particle and debited
 *   from the other. It lives in the force loop, not the density loop, for
 *   two reasons that are one: the force loop's dispatch fires BOTH sides
 *   whenever EITHER kernel reaches, so the mirrored credit/debit pair is
 *   never split (a density-loop placement fabricates energy at
 *   h_i != h_j); and the loop runs after the extra ghost has set this
 *   step's coefficient and before cooling reads `u`, so the trigger acts
 *   within the step it fires. SPHENIX's own artificial viscosity lives in
 *   the force loop for the identical reason. Closure-independent.
 *
 * All three operators need a stable per-particle density: the density loop
 * here runs interleaved with SPH's own density accumulation, so `p->rho` is
 * a partial sum, not a density, at the point those pairwise calls run. All
 * therefore read `p->feedback_data.rho_prev`, a comoving density snapshot
 * cached once per step by `radiation_snapshot_part_propagation` (before the
 * per-step density-accumulator reset), the same snapshot for every loop.
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
 * @param div_F_i (return, accumulated) Particle i's div(F) accumulator.
 * @param div_F_j (return, accumulated) Particle j's div(F) accumulator.
 */
__attribute__((always_inline)) INLINE static void
radiation_divergence_accumulate_band(const float dx[3], float r_inv,
                                     float wi_dr, float wj_dr, float mi,
                                     float mj, float rho_i, float rho_j,
                                     const float F_i[3], const float F_j[3],
                                     float *div_F_i, float *div_F_j) {

  const float Fi_dot_dx = F_i[0] * dx[0] + F_i[1] * dx[1] + F_i[2] * dx[2];
  const float Fj_dot_dx = F_j[0] * dx[0] + F_j[1] * dx[1] + F_j[2] * dx[2];

  const float Phi_ij =
      Fi_dot_dx / rho_i * wi_dr * r_inv + Fj_dot_dx / rho_j * wj_dr * r_inv;

  *div_F_i += mj * Phi_ij;
  *div_F_j += -mi * Phi_ij;
}

/**
 * @brief Stage-2 van Leer slope limiter for the midpoint reconstruction of
 * the Stage-1 jump (design-lw-fuv-design-b-dissipation.md Section 5.2), in
 * MAGMA2's compiled form (`origin/MAGMA2_matthieu:src/hydro/MAGMA/
 * hydro_iact.h`, Rosswog 2020b Eq. 21-23) transcribed to `u_V = rho*u`.
 *
 * Returns 0 for opposite-sign gradients, so a single-particle dip (the
 * configuration Stage 1 exists for) keeps the raw jump. The `A_ij = -1`
 * pole is returned as 0 rather than MAGMA's own 1: an optimized build
 * carries `-ffast-math`, so the clamp cannot be relied on to turn the
 * resulting infinity back into a finite limiter value.
 *
 * @param dx Comoving separation vector (pi - pj).
 * @param r Comoving particle separation.
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param g_i Particle i's previous-step `grad(u_V)` (this band).
 * @param g_j Particle j's previous-step `grad(u_V)` (this band).
 * @return The limiter `Phi_ij`, in [0, 1].
 */
__attribute__((always_inline)) INLINE static float
radiation_dissipation_van_leer_limiter(const float dx[3], float r, float hi,
                                       float hj, const float g_i[3],
                                       const float g_j[3]) {

  const float A_num = g_i[0] * dx[0] + g_i[1] * dx[1] + g_i[2] * dx[2];
  const float A_den = g_j[0] * dx[0] + g_j[1] * dx[1] + g_j[2] * dx[2];
  const float A_ij = (A_den != 0.f) ? A_num / A_den : 0.f;

  const float one_plus_A = 1.f + A_ij;
  const float A_denominator = one_plus_A * one_plus_A;
  const float fraction =
      (A_denominator > 0.f) ? 4.f * A_ij / A_denominator : 0.f;
  const float fraction_capped = min(fraction, 1.f);
  const float limiter = max(fraction_capped, 0.f);

  const float eta_ij = r / max(hi, hj);
  const float d_eta = eta_ij - RADIATION_LW_FUV_DISSIPATION_ETA_CRIT;
  const float exp_term = (eta_ij < RADIATION_LW_FUV_DISSIPATION_ETA_CRIT)
                             ? expf(-25.f * d_eta * d_eta)
                             : 1.f;

  return limiter * exp_term;
}

/**
 * @brief Band-specific pairwise contribution to each particle's kernel-mean
 * `|rho_prev*u_prev|` reference accumulator, the local field scale the
 * Stage-1 negativity trigger divides an undershoot by (radiation_isrf.c's
 * #radiation_update_dissipation_alpha_band).
 *
 * Accumulated in the density loop, from the stable `u_*_prev` snapshot and
 * `rho_prev`, so the value the trigger reads does not drift across a
 * particle's h-iterations. The dissipation term the trigger drives is
 * accumulated separately, in the force loop
 * (#radiation_dissipation_force_accumulate_band).
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
 * @brief Band-specific pairwise contribution to particle i's Stage-1
 * artificial-dissipation source term (design-lw-fuv-design-b-
 * dissipation.md Section 3.1), and mirrored (mass-weighted, opposite
 * sign) contribution to particle j's.
 *
 * `v_sig,ij = max(alpha_i, alpha_j) * min(c_hyp_i, c_hyp_j)` is a signal
 * VELOCITY, with no `h` factor: the length scale enters only through
 * `Wbar_ij`'s own `h^-(dim+1)` normalisation, exactly as for an ordinary
 * SPH Laplacian.
 *
 * Accumulated in the force loop, which runs after the density ghost has
 * produced the intermediate state `u* = e*u_prev + dt*phi*(source - div_F)`
 * and after the extra ghost has set this step's `alpha`: the jump is
 * therefore built from the LIVE `u_FUV`/`u_LW` (`u*`), not from the
 * `u_*_prev` snapshot the density loop needs. The force loop runs exactly
 * once per step, so there is no h-iteration stability requirement here.
 *
 * The credit/debit pair is applied unconditionally, with no mutual-reach
 * gate: the force loop's dispatch fires both sides whenever either kernel
 * reaches, so the pair is always mirrored and `sum_i m_i*dissipation_u_i`
 * is exactly zero over a same-bin active-active pair for any h_i, h_j.
 * `Psi_ji = -Psi_ij` under the exchange `i <-> j`, `dx -> -dx`: `Wbar_ij`
 * is symmetric, `d_ij` is antisymmetric, and (Stage 2) the van Leer
 * limiter is invariant under `A_ij -> 1/A_ij`, so the reconstruction does
 * not break the antisymmetry either.
 *
 * Under #RADIATION_LW_FUV_DISSIPATION_RECONSTRUCTION (Stage 2) the jump is
 * first reconstructed to the pair midpoint with the limiter above, using
 * each particle's #grad_u_FUV_prev/LW_prev. That field is written at the
 * end of the extra ghost, which precedes this loop, so once accumulated
 * here the reconstruction reads THIS step's finalized gradient rather than
 * the previous step's.
 *
 * @param dx Comoving separation vector (pi - pj).
 * @param r Comoving particle separation.
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param wi_dr See #radiation_divergence_accumulate_band.
 * @param wj_dr See #radiation_divergence_accumulate_band.
 * @param mi Particle i's mass.
 * @param mj Particle j's mass.
 * @param rho_i Particle i's cached comoving density snapshot.
 * @param rho_j Particle j's cached comoving density snapshot.
 * @param c_hyp_i Particle i's own hyperbolic propagation speed.
 * @param c_hyp_j Particle j's own hyperbolic propagation speed.
 * @param alpha_i Particle i's own dissipation coefficient (this band).
 * @param alpha_j Particle j's own dissipation coefficient (this band).
 * @param u_i Particle i's live specific field `u*` (this band).
 * @param u_j Particle j's live specific field `u*` (this band).
 * @param grad_u_i_prev Particle i's #grad_u_FUV_prev/LW_prev (this band).
 * @param grad_u_j_prev Particle j's #grad_u_FUV_prev/LW_prev (this band).
 * @param dissipation_u_i (return, accumulated) Particle i's dissipation
 * source-term accumulator.
 * @param dissipation_u_j (return, accumulated) Particle j's dissipation
 * source-term accumulator.
 */
__attribute__((always_inline)) INLINE static void
radiation_dissipation_force_accumulate_band(
    const float dx[3], float r, float hi, float hj, float wi_dr, float wj_dr,
    float mi, float mj, float rho_i, float rho_j, float c_hyp_i, float c_hyp_j,
    float alpha_i, float alpha_j, float u_i, float u_j,
    const float grad_u_i_prev[3], const float grad_u_j_prev[3],
    float *dissipation_u_i, float *dissipation_u_j) {

#ifdef RADIATION_LW_FUV_DISSIPATION_RECONSTRUCTION
  const int use_reconstruction = 1;
#else
  const int use_reconstruction = 0;
#endif

  float d_ij = rho_i * u_i - rho_j * u_j;

  if (use_reconstruction) {
    const float g_i[3] = {rho_i * grad_u_i_prev[0], rho_i * grad_u_i_prev[1],
                          rho_i * grad_u_i_prev[2]};
    const float g_j[3] = {rho_j * grad_u_j_prev[0], rho_j * grad_u_j_prev[1],
                          rho_j * grad_u_j_prev[2]};
    const float Phi_ij =
        radiation_dissipation_van_leer_limiter(dx, r, hi, hj, g_i, g_j);
    const float g_sum_dot_dx = (g_i[0] + g_j[0]) * dx[0] +
                               (g_i[1] + g_j[1]) * dx[1] +
                               (g_i[2] + g_j[2]) * dx[2];
    d_ij -= Phi_ij * 0.5f * g_sum_dot_dx;
  }

  const float Wbar_ij = 0.5f * (wi_dr + wj_dr);
  const float v_sig_ij = max(alpha_i, alpha_j) * min(c_hyp_i, c_hyp_j);
  const float Psi_ij = v_sig_ij * d_ij * Wbar_ij / (rho_i * rho_j);

  *dissipation_u_i += mj * Psi_ij;
  *dissipation_u_j += -mi * Psi_ij;
}

/**
 * @brief M1 closure tensor `D(f)` for one particle, one band, built from its
 * own `(u, F, c_M)` (design-lw-fuv-m1-upgrade.md "New pieces"). `c_M` is the
 * same speed already carried as #feedback_part_data.c_hyp (D3: reinterpreted
 * as the fastest M1 characteristic, `f=1`, not a new field).
 *
 * `f = min(1, |F|/(c_M*u))` for `u > 0`, `f = 0` for `u <= 0`;
 * `chi(f) = (3+4f^2)/(5+2*sqrt(4-3f^2))`;
 * `D(f) = (1-chi)/2 I + (3chi-1)/2 (n dyadic n)`, `n = F/|F|`.
 *
 * Zero-flux guard, mandatory: `F = 0` is every particle's initial condition
 * and permanent far-field state, not a corner case. `F2 = F.F`,
 * `F_inv = (F2 > 0) ? 1/sqrt(F2) : 0`, `n = F*F_inv` -- the same convention
 * #radiation_flux_dissipation_accumulate_band already uses for `F_inv_i`.
 * `f`'s own division is guarded the same way: `c_M*u` is computed once and
 * only divided into when it is strictly positive, which also folds in the
 * `u <= 0` case (`f = 0`) without a separate branch. At `F = 0`, `f = 0`,
 * `chi = 1/3`, the `(3*chi-1)/2 = 0` coefficient multiplies the guarded,
 * well-defined zero `n` rather than a NaN.
 *
 * @param u This band's ghost-finalized specific field for this particle.
 * @param F This particle's tracked flux (this band).
 * @param c_M This particle's own #feedback_part_data.c_hyp.
 * @param D (return) The 3x3 closure tensor.
 */
__attribute__((always_inline)) INLINE static void
radiation_get_m1_closure_tensor_band(float u, const float F[3], float c_M,
                                     float D[3][3]) {

  const float F2 = F[0] * F[0] + F[1] * F[1] + F[2] * F[2];
  const float F_inv = (F2 > 0.f) ? 1.f / sqrtf(F2) : 0.f;
  const float Fmag = F2 * F_inv; /* sqrt(F2), no second sqrtf call */
  const float n[3] = {F[0] * F_inv, F[1] * F_inv, F[2] * F_inv};

  const float denom = c_M * u;
  const float f = (denom > 0.f) ? min(Fmag / denom, 1.f) : 0.f;

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
 * `diffmode == 2` form specifically (design-lw-fuv-m1-upgrade.md, pinned by
 * plan review): `tempi - tempj` on a shared, averaged kernel derivative
 * `(wi_dr + wj_dr)*0.5`, matching `src/rt/SPHM1RT/rt_gradients.h`'s
 * `radiation_gradient_aniso_SPH` `diffmode==2` branch
 * (`src/rt/SPHM1RT/rt_iact.h:582-627` calls it with `diffmodeaniso = 2`) and
 * this project's own existing `d_ij = rho_i*u_i - rho_j*u_j` jump
 * construction. `D_i`, `D_j` reduce to `(1/3) I` at `f=0` (both particles'
 * fluxes zero, the isotropic P1 limit), so this reduces to the old scalar
 * form's structure with the `1/3` now explicit rather than folded away.
 *
 * @param dx Comoving separation vector (pi - pj).
 * @param r_inv Inverse comoving particle separation.
 * @param wi_dr See #radiation_divergence_accumulate_band.
 * @param wj_dr See #radiation_divergence_accumulate_band.
 * @param mi Particle i's mass.
 * @param mj Particle j's mass.
 * @param rho_i Particle i's cached comoving density snapshot.
 * @param rho_j Particle j's cached comoving density snapshot.
 * @param u_i Particle i's ghost-finalized specific field (this band).
 * @param u_j Particle j's ghost-finalized specific field (this band).
 * @param D_i Particle i's own M1 closure tensor (this band), from
 * #radiation_get_m1_closure_tensor_band.
 * @param D_j Particle j's own M1 closure tensor (this band).
 * @param grad_u_i (return, accumulated) Particle i's grad(u) accumulator.
 * @param grad_u_j (return, accumulated) Particle j's grad(u) accumulator.
 */
__attribute__((always_inline)) INLINE static void
radiation_gradient_accumulate_band(const float dx[3], float r_inv, float wi_dr,
                                   float wj_dr, float mi, float mj, float rho_i,
                                   float rho_j, float u_i, float u_j,
                                   const float D_i[3][3], const float D_j[3][3],
                                   float grad_u_i[3], float grad_u_j[3]) {

  const float rho_i_inv = 1.f / rho_i;
  const float rho_j_inv = 1.f / rho_j;
  const float wbar_dr = 0.5f * (wi_dr + wj_dr);

  float temp_i[3], temp_j[3];
  for (int k = 0; k < 3; k++) {
    const float Di_dot_dx =
        D_i[k][0] * dx[0] + D_i[k][1] * dx[1] + D_i[k][2] * dx[2];
    const float Dj_dot_dx =
        D_j[k][0] * dx[0] + D_j[k][1] * dx[1] + D_j[k][2] * dx[2];
    temp_i[k] = Di_dot_dx * rho_i * u_i * r_inv;
    temp_j[k] = Dj_dot_dx * rho_j * u_j * r_inv;
  }

  const float fac_i = mj * rho_i_inv * rho_i_inv * wbar_dr;
  const float fac_j = mi * rho_j_inv * rho_j_inv * wbar_dr;

  for (int k = 0; k < 3; k++) {
    grad_u_i[k] += -(temp_i[k] - temp_j[k]) * fac_i;
    grad_u_j[k] += -(temp_i[k] - temp_j[k]) * fac_j;
  }
}

/**
 * @brief Band-specific pairwise contribution to the Stage-3 anisotropic
 * flux-dissipation source term (design-lw-fuv-design-b-dissipation.md
 * Section 5.2), the `diffmode==2` form of
 * `src/rt/SPHM1RT/rt_gradients.h`'s `radiation_gradient_aniso_SPH` applied
 * to `D_f * psi` along `n n`.
 *
 * Accumulated in the gradient loop, which already runs after
 * `psi = div(F)` is final in the density ghost, so no third loop is needed.
 * NOT antisymmetric: the flux is not a conserved sum, so both particles
 * take the same-signed shared scalar and each divides by its own density.
 * A particle whose flux is exactly zero (every particle's initial state,
 * and the permanent far-field state) has no direction to be anisotropic
 * along and contributes nothing.
 *
 * The pair's signal velocity is `min(c_hyp_i, c_hyp_j)` as for the Stage-1
 * term (Section 3.2), while `alpha_f` and `h` stay per particle.
 *
 * @param dx Comoving separation vector (pi - pj).
 * @param r_inv Inverse comoving particle separation.
 * @param wi_dr See #radiation_divergence_accumulate_band.
 * @param wj_dr See #radiation_divergence_accumulate_band.
 * @param mi Particle i's mass.
 * @param mj Particle j's mass.
 * @param rho_i Particle i's cached comoving density snapshot.
 * @param rho_j Particle j's cached comoving density snapshot.
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param c_hyp_i Particle i's own hyperbolic propagation speed.
 * @param c_hyp_j Particle j's own hyperbolic propagation speed.
 * @param alpha_f_i Particle i's own flux-dissipation coefficient (band).
 * @param alpha_f_j Particle j's own flux-dissipation coefficient (band).
 * @param F_i Particle i's tracked flux (this band).
 * @param F_j Particle j's tracked flux (this band).
 * @param psi_i Particle i's finalized `div(F)` (this band).
 * @param psi_j Particle j's finalized `div(F)` (this band).
 * @param dissipation_F_i (return, accumulated) Particle i's flux-dissipation
 * source-term accumulator.
 * @param dissipation_F_j (return, accumulated) Particle j's flux-dissipation
 * source-term accumulator.
 */
__attribute__((always_inline)) INLINE static void
radiation_flux_dissipation_accumulate_band(
    const float dx[3], float r_inv, float wi_dr, float wj_dr, float mi,
    float mj, float rho_i, float rho_j, float hi, float hj, float c_hyp_i,
    float c_hyp_j, float alpha_f_i, float alpha_f_j, const float F_i[3],
    const float F_j[3], float psi_i, float psi_j, float dissipation_F_i[3],
    float dissipation_F_j[3]) {

  const float v_sig = min(c_hyp_i, c_hyp_j);

  const float F2_i = F_i[0] * F_i[0] + F_i[1] * F_i[1] + F_i[2] * F_i[2];
  const float F2_j = F_j[0] * F_j[0] + F_j[1] * F_j[1] + F_j[2] * F_j[2];
  const float F_inv_i = (F2_i > 0.f) ? 1.f / sqrtf(F2_i) : 0.f;
  const float F_inv_j = (F2_j > 0.f) ? 1.f / sqrtf(F2_j) : 0.f;

  /* `n_i . dx`, and the scalar that turns F_i[k] into t_i[k] once
     multiplied by it: t_i = rho_i * D_f,i * psi_i * n_i * (n_i . dx). */
  const float ni_dot_dx =
      (F_i[0] * dx[0] + F_i[1] * dx[1] + F_i[2] * dx[2]) * F_inv_i;
  const float nj_dot_dx =
      (F_j[0] * dx[0] + F_j[1] * dx[1] + F_j[2] * dx[2]) * F_inv_j;
  const float t_fac_i =
      rho_i * alpha_f_i * v_sig * hi * psi_i * ni_dot_dx * F_inv_i;
  const float t_fac_j =
      rho_j * alpha_f_j * v_sig * hj * psi_j * nj_dot_dx * F_inv_j;

  const float Wbar_ij = 0.5f * (wi_dr + wj_dr);
  const float weight_i = -mj * Wbar_ij * r_inv / (rho_i * rho_i);
  const float weight_j = -mi * Wbar_ij * r_inv / (rho_j * rho_j);

  for (int k = 0; k < 3; k++) {
    const float t_diff = t_fac_i * F_i[k] - t_fac_j * F_j[k];
    dissipation_F_i[k] += weight_i * t_diff;
    dissipation_F_j[k] += weight_j * t_diff;
  }
}

/**
 * @brief `div(F)` propagation interaction between two particles
 * (symmetric): both particles' accumulators are updated.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param us Unit system (unused: the SPH operator needs only positions,
 * masses, the cached density snapshot, and the tracked flux).
 */
__attribute__((always_inline)) INLINE static void runner_iact_isrf_propagation(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const struct unit_system *us) {

  const float r = sqrtf(r2);
  const float r_inv = r ? 1.f / r : 0.f;
  const float hi_inv = 1.f / hi;
  const float hj_inv = 1.f / hj;

  float wi, wi_dx, wj, wj_dx;
  kernel_deval(r * hi_inv, &wi, &wi_dx);
  kernel_deval(r * hj_inv, &wj, &wj_dx);
  const float wi_dr = wi_dx * pow_dimension_plus_one(hi_inv);
  const float wj_dr = wj_dx * pow_dimension_plus_one(hj_inv);
  wi *= pow_dimension(hi_inv);
  wj *= pow_dimension(hj_inv);

  struct feedback_part_data *fdi = &pi->feedback_data;
  struct feedback_part_data *fdj = &pj->feedback_data;
  const float rho_i = fdi->rho_prev;
  const float rho_j = fdj->rho_prev;
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  radiation_divergence_accumulate_band(
      dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, fdi->specific_flux_FUV,
      fdj->specific_flux_FUV, &fdi->div_specific_flux_FUV,
      &fdj->div_specific_flux_FUV);
  radiation_divergence_accumulate_band(
      dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, fdi->specific_flux_LW,
      fdj->specific_flux_LW, &fdi->div_specific_flux_LW,
      &fdj->div_specific_flux_LW);

  radiation_dissipation_reference_accumulate_band(
      wi, wj, mi, mj, rho_i, rho_j, fdi->u_FUV_prev, fdj->u_FUV_prev,
      &fdi->ngb_mean_abs_u_V_FUV, &fdj->ngb_mean_abs_u_V_FUV);
  radiation_dissipation_reference_accumulate_band(
      wi, wj, mi, mj, rho_i, rho_j, fdi->u_LW_prev, fdj->u_LW_prev,
      &fdi->ngb_mean_abs_u_V_LW, &fdj->ngb_mean_abs_u_V_LW);
}

/**
 * @brief `div(F)` propagation interaction between two particles
 * (non-symmetric): only particle i's accumulators are updated.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (its own accumulators not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
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
  const float r_inv = r ? 1.f / r : 0.f;
  const float hi_inv = 1.f / hi;
  const float hj_inv = 1.f / hj;

  float wi, wi_dx, wj, wj_dx;
  kernel_deval(r * hi_inv, &wi, &wi_dx);
  kernel_deval(r * hj_inv, &wj, &wj_dx);
  const float wi_dr = wi_dx * pow_dimension_plus_one(hi_inv);
  const float wj_dr = wj_dx * pow_dimension_plus_one(hj_inv);
  wi *= pow_dimension(hi_inv);
  wj *= pow_dimension(hj_inv);

  struct feedback_part_data *fdi = &pi->feedback_data;
  const struct feedback_part_data *fdj = &pj->feedback_data;
  const float rho_i = fdi->rho_prev;
  const float rho_j = fdj->rho_prev;
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  /* Particle j's own accumulator is not touched (non-symmetric): pass a
   * discarded local, seeded to 0 rather than read from fdj, as the
   * required (return, accumulated) output. */
  float unused_div_specific_flux_FUV = 0.f;
  float unused_div_specific_flux_LW = 0.f;
  float unused_ngb_mean_abs_u_V_FUV = 0.f;
  float unused_ngb_mean_abs_u_V_LW = 0.f;

  radiation_divergence_accumulate_band(
      dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, fdi->specific_flux_FUV,
      fdj->specific_flux_FUV, &fdi->div_specific_flux_FUV,
      &unused_div_specific_flux_FUV);
  radiation_divergence_accumulate_band(
      dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, fdi->specific_flux_LW,
      fdj->specific_flux_LW, &fdi->div_specific_flux_LW,
      &unused_div_specific_flux_LW);

  radiation_dissipation_reference_accumulate_band(
      wi, wj, mi, mj, rho_i, rho_j, fdi->u_FUV_prev, fdj->u_FUV_prev,
      &fdi->ngb_mean_abs_u_V_FUV, &unused_ngb_mean_abs_u_V_FUV);
  radiation_dissipation_reference_accumulate_band(
      wi, wj, mi, mj, rho_i, rho_j, fdi->u_LW_prev, fdj->u_LW_prev,
      &fdi->ngb_mean_abs_u_V_LW, &unused_ngb_mean_abs_u_V_LW);
}

/**
 * @brief `grad(u)` interaction between two particles (symmetric): both
 * particles' accumulators are updated.
 *
 * Runs in the gradient loop, after the density ghost has finalized `u_FUV`/
 * `u_LW` for this step (the exact-relaxation `u` update, radiation_isrf.c):
 * reads them directly, not a `_prev` snapshot. `u_FUV`/`u_LW` are written
 * again later in this same step, by the end-force ghost's Stage-1
 * dissipation correction (radiation_isrf.c's
 * #radiation_end_force_propagation), before star feedback injection ever
 * runs.
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

  float D_FUV_i[3][3], D_FUV_j[3][3], D_LW_i[3][3], D_LW_j[3][3];
  radiation_get_m1_closure_tensor_band(fdi->u_FUV, fdi->specific_flux_FUV,
                                       fdi->c_hyp, D_FUV_i);
  radiation_get_m1_closure_tensor_band(fdj->u_FUV, fdj->specific_flux_FUV,
                                       fdj->c_hyp, D_FUV_j);
  radiation_get_m1_closure_tensor_band(fdi->u_LW, fdi->specific_flux_LW,
                                       fdi->c_hyp, D_LW_i);
  radiation_get_m1_closure_tensor_band(fdj->u_LW, fdj->specific_flux_LW,
                                       fdj->c_hyp, D_LW_j);

  radiation_gradient_accumulate_band(dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i,
                                     rho_j, fdi->u_FUV, fdj->u_FUV, D_FUV_i,
                                     D_FUV_j, fdi->grad_u_FUV, fdj->grad_u_FUV);
  radiation_gradient_accumulate_band(dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i,
                                     rho_j, fdi->u_LW, fdj->u_LW, D_LW_i,
                                     D_LW_j, fdi->grad_u_LW, fdj->grad_u_LW);

  /* Stage 3 owns per-particle fields that only exist when the stage is
     built, so its call sites are guarded rather than gated on a runtime
     flag; the formulas above stay compiled in either state. */
#ifdef RADIATION_LW_FUV_DISSIPATION_ANISOTROPIC_FLUX
  radiation_flux_dissipation_accumulate_band(
      dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, hi, hj, fdi->c_hyp,
      fdj->c_hyp, fdi->dissipation_alpha_flux_FUV,
      fdj->dissipation_alpha_flux_FUV, fdi->specific_flux_FUV,
      fdj->specific_flux_FUV, fdi->div_specific_flux_FUV,
      fdj->div_specific_flux_FUV, fdi->dissipation_F_FUV,
      fdj->dissipation_F_FUV);
  radiation_flux_dissipation_accumulate_band(
      dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, hi, hj, fdi->c_hyp,
      fdj->c_hyp, fdi->dissipation_alpha_flux_LW,
      fdj->dissipation_alpha_flux_LW, fdi->specific_flux_LW,
      fdj->specific_flux_LW, fdi->div_specific_flux_LW,
      fdj->div_specific_flux_LW, fdi->dissipation_F_LW, fdj->dissipation_F_LW);
#endif
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

  /* Particle j is `const` here (non-symmetric): its own accumulator is not
   * touched, so pass a discarded, zero-seeded local as the writable
   * destination the shared accumulator function requires for j. */
  float unused_grad_u_FUV[3] = {0.f, 0.f, 0.f};
  float unused_grad_u_LW[3] = {0.f, 0.f, 0.f};

  float D_FUV_i[3][3], D_FUV_j[3][3], D_LW_i[3][3], D_LW_j[3][3];
  radiation_get_m1_closure_tensor_band(fdi->u_FUV, fdi->specific_flux_FUV,
                                       fdi->c_hyp, D_FUV_i);
  radiation_get_m1_closure_tensor_band(fdj->u_FUV, fdj->specific_flux_FUV,
                                       fdj->c_hyp, D_FUV_j);
  radiation_get_m1_closure_tensor_band(fdi->u_LW, fdi->specific_flux_LW,
                                       fdi->c_hyp, D_LW_i);
  radiation_get_m1_closure_tensor_band(fdj->u_LW, fdj->specific_flux_LW,
                                       fdj->c_hyp, D_LW_j);

  radiation_gradient_accumulate_band(
      dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, fdi->u_FUV, fdj->u_FUV,
      D_FUV_i, D_FUV_j, fdi->grad_u_FUV, unused_grad_u_FUV);
  radiation_gradient_accumulate_band(dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i,
                                     rho_j, fdi->u_LW, fdj->u_LW, D_LW_i,
                                     D_LW_j, fdi->grad_u_LW, unused_grad_u_LW);

  /* See the symmetric variant above for why this is guarded rather than
     gated on a runtime flag. */
#ifdef RADIATION_LW_FUV_DISSIPATION_ANISOTROPIC_FLUX
  float unused_dissipation_F_FUV[3] = {0.f, 0.f, 0.f};
  float unused_dissipation_F_LW[3] = {0.f, 0.f, 0.f};

  radiation_flux_dissipation_accumulate_band(
      dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, hi, hj, fdi->c_hyp,
      fdj->c_hyp, fdi->dissipation_alpha_flux_FUV,
      fdj->dissipation_alpha_flux_FUV, fdi->specific_flux_FUV,
      fdj->specific_flux_FUV, fdi->div_specific_flux_FUV,
      fdj->div_specific_flux_FUV, fdi->dissipation_F_FUV,
      unused_dissipation_F_FUV);
  radiation_flux_dissipation_accumulate_band(
      dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, hi, hj, fdi->c_hyp,
      fdj->c_hyp, fdi->dissipation_alpha_flux_LW,
      fdj->dissipation_alpha_flux_LW, fdi->specific_flux_LW,
      fdj->specific_flux_LW, fdi->div_specific_flux_LW,
      fdj->div_specific_flux_LW, fdi->dissipation_F_LW,
      unused_dissipation_F_LW);
#endif
}

/**
 * @brief Stage-1 artificial-dissipation interaction between two particles
 * (symmetric): both particles' accumulators are updated.
 *
 * Runs in the force loop, after the density ghost has produced `u*` and the
 * extra ghost has set this step's #dissipation_alpha_FUV/LW: reads the live
 * `u_FUV`/`u_LW` directly, not a `_prev` snapshot. The force loop's
 * dispatch fires both sides of a pair whenever either kernel reaches, which
 * is what keeps the mirrored credit/debit pair whole at h_i != h_j; see
 * #radiation_dissipation_force_accumulate_band.
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

  /* The Stage-2 gradients only exist when the stage is built, so the two
     pointers are selected here; the formula they feed stays compiled in
     both states. */
#ifdef RADIATION_LW_FUV_DISSIPATION_RECONSTRUCTION
  const float *const g_FUV_i = fdi->grad_u_FUV_prev;
  const float *const g_FUV_j = fdj->grad_u_FUV_prev;
  const float *const g_LW_i = fdi->grad_u_LW_prev;
  const float *const g_LW_j = fdj->grad_u_LW_prev;
#else
  const float g_absent[3] = {0.f, 0.f, 0.f};
  const float *const g_FUV_i = g_absent;
  const float *const g_FUV_j = g_absent;
  const float *const g_LW_i = g_absent;
  const float *const g_LW_j = g_absent;
#endif

  radiation_dissipation_force_accumulate_band(
      dx, r, hi, hj, wi_dr, wj_dr, mi, mj, rho_i, rho_j, fdi->c_hyp, fdj->c_hyp,
      fdi->dissipation_alpha_FUV, fdj->dissipation_alpha_FUV, fdi->u_FUV,
      fdj->u_FUV, g_FUV_i, g_FUV_j, &fdi->dissipation_u_FUV,
      &fdj->dissipation_u_FUV);
  radiation_dissipation_force_accumulate_band(
      dx, r, hi, hj, wi_dr, wj_dr, mi, mj, rho_i, rho_j, fdi->c_hyp, fdj->c_hyp,
      fdi->dissipation_alpha_LW, fdj->dissipation_alpha_LW, fdi->u_LW,
      fdj->u_LW, g_LW_i, g_LW_j, &fdi->dissipation_u_LW,
      &fdj->dissipation_u_LW);
}

/**
 * @brief Stage-1 artificial-dissipation interaction between two particles
 * (non-symmetric): only particle i's accumulator is updated.
 *
 * The force loop reaches this variant once per side, so a pair whose two
 * sides are both dispatched still receives the mirrored credit/debit pair
 * in full; see #runner_iact_isrf_dissipation.
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

  /* Particle j's own accumulator is not touched (non-symmetric): pass a
   * discarded local, seeded to 0 rather than read from fdj, as the
   * required (return, accumulated) output. */
  float unused_dissipation_u_FUV = 0.f;
  float unused_dissipation_u_LW = 0.f;

  /* See the symmetric variant above for why these pointers are selected
     here rather than guarding the call itself. */
#ifdef RADIATION_LW_FUV_DISSIPATION_RECONSTRUCTION
  const float *const g_FUV_i = fdi->grad_u_FUV_prev;
  const float *const g_FUV_j = fdj->grad_u_FUV_prev;
  const float *const g_LW_i = fdi->grad_u_LW_prev;
  const float *const g_LW_j = fdj->grad_u_LW_prev;
#else
  const float g_absent[3] = {0.f, 0.f, 0.f};
  const float *const g_FUV_i = g_absent;
  const float *const g_FUV_j = g_absent;
  const float *const g_LW_i = g_absent;
  const float *const g_LW_j = g_absent;
#endif

  radiation_dissipation_force_accumulate_band(
      dx, r, hi, hj, wi_dr, wj_dr, mi, mj, rho_i, rho_j, fdi->c_hyp, fdj->c_hyp,
      fdi->dissipation_alpha_FUV, fdj->dissipation_alpha_FUV, fdi->u_FUV,
      fdj->u_FUV, g_FUV_i, g_FUV_j, &fdi->dissipation_u_FUV,
      &unused_dissipation_u_FUV);
  radiation_dissipation_force_accumulate_band(
      dx, r, hi, hj, wi_dr, wj_dr, mi, mj, rho_i, rho_j, fdi->c_hyp, fdj->c_hyp,
      fdi->dissipation_alpha_LW, fdj->dissipation_alpha_LW, fdi->u_LW,
      fdj->u_LW, g_LW_i, g_LW_j, &fdi->dissipation_u_LW,
      &unused_dissipation_u_LW);
}

#endif /* SWIFT_RADIATION_PROPAGATION_IACT_GEAR_H */
