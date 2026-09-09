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
 * @brief Gas-gas density-loop and gradient-loop hooks for the hyperbolic
 * P1-relaxation propagation of the u_FUV/specific_flux_FUV (and
 * u_LW/specific_flux_LW) fields.
 *
 * Two pairwise SPH operators are accumulated here, in two different loops:
 *
 * - `div(F)` (density loop, this file's
 * `runner_iact_[nonsym_]isrf_propagation`): the shared-coefficient construction
 * mirroring `src/rt/SPHM1RT/rt_gradients.h`'s `radiation_divergence_SPH`
 *   `diffmode==1` branch. Exactly mass-conserving under transport alone for
 *   any h_i != h_j, rho_i != rho_j (a single shared scalar built from both
 *   particles' own kernel-gradient terms, applied with mirrored mass/sign to
 *   each side).
 * - `grad(u)` (gradient loop, `runner_iact_[nonsym_]isrf_gradient`): the
 *   difference-on-(rho*u) construction mirroring the same file's
 *   `radiation_gradient_SPH` `diffmode==0` branch. This is the operator that
 *   is minus the adjoint of the diffmode==1 divergence above in the m*rho
 *   inner product, which the staggered exact-relaxation time integrator
 *   (radiation_isrf.c) needs for stability on a disordered particle
 *   distribution.
 *
 * Both operators need a stable per-particle density: the density loop here
 * runs interleaved with SPH's own density accumulation, so `p->rho` is a
 * partial sum, not a density, at the point these pairwise calls run. Both
 * therefore read `p->feedback_data.rho_prev`, a comoving density snapshot
 * cached once per step by `radiation_snapshot_part_propagation` (before the
 * per-step density-accumulator reset), the same snapshot for both loops so
 * `div` and `grad` are built from the identical `rho_i`, `rho_j` values
 * their skew-adjoint pairing requires.
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
 * @brief Band-specific pairwise contribution to particle i's Stage-1
 * artificial-dissipation source term (design-lw-fuv-design-b-
 * dissipation.md Section 3.1), and mirrored (mass-weighted, opposite
 * sign) contribution to particle j's, plus each particle's own kernel-mean
 * `|rho_prev*u_prev|` reference accumulator the negativity trigger divides
 * by (radiation_isrf.c's #radiation_update_dissipation_alpha_band).
 *
 * `v_sig,ij = max(alpha_i, alpha_j) * min(c_hyp_i, c_hyp_j)` is a signal
 * VELOCITY, with no `h` factor: the length scale enters only through
 * `Wbar_ij`'s own `h^-(dim+1)` normalisation, exactly as for an ordinary
 * SPH Laplacian. `u_i_prev`/`u_j_prev` (not the live, in-progress `u`) and
 * `rho_i`/`rho_j` (#feedback_part_data.rho_prev) keep this term stable
 * across a particle's h-iterations, as #radiation_gradient_accumulate_band
 * below already requires for `d_ij`.
 *
 * Under #RADIATION_LW_FUV_DISSIPATION_RECONSTRUCTION (Stage 2) the jump is
 * first reconstructed to the pair midpoint with the limiter above, using
 * each particle's previous-step gradient.
 *
 * @param dx Comoving separation vector (pi - pj).
 * @param r Comoving particle separation.
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param wi Particle i's own kernel value, W(r/h_i)*h_i^-dim.
 * @param wj Particle j's own kernel value, W(r/h_j)*h_j^-dim.
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
 * @param u_i_prev Particle i's snapshotted specific field (this band).
 * @param u_j_prev Particle j's snapshotted specific field (this band).
 * @param grad_u_i_prev Particle i's previous-step `grad(u)` (this band).
 * @param grad_u_j_prev Particle j's previous-step `grad(u)` (this band).
 * @param dissipation_u_i (return, accumulated) Particle i's dissipation
 * source-term accumulator.
 * @param dissipation_u_j (return, accumulated) Particle j's dissipation
 * source-term accumulator.
 * @param ngb_mean_abs_u_V_i (return, accumulated) Particle i's kernel-mean
 * `|rho_prev*u_prev|` accumulator.
 * @param ngb_mean_abs_u_V_j (return, accumulated) Particle j's kernel-mean
 * `|rho_prev*u_prev|` accumulator.
 */
__attribute__((always_inline)) INLINE static void
radiation_dissipation_accumulate_band(
    const float dx[3], float r, float hi, float hj, float wi, float wj,
    float wi_dr, float wj_dr, float mi, float mj, float rho_i, float rho_j,
    float c_hyp_i, float c_hyp_j, float alpha_i, float alpha_j, float u_i_prev,
    float u_j_prev, const float grad_u_i_prev[3], const float grad_u_j_prev[3],
    float *dissipation_u_i, float *dissipation_u_j, float *ngb_mean_abs_u_V_i,
    float *ngb_mean_abs_u_V_j) {

#ifdef RADIATION_LW_FUV_DISSIPATION_RECONSTRUCTION
  const int use_reconstruction = 1;
#else
  const int use_reconstruction = 0;
#endif

  float d_ij = rho_i * u_i_prev - rho_j * u_j_prev;

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

  *ngb_mean_abs_u_V_i += (mj / rho_j) * wi * fabsf(rho_j * u_j_prev);
  *ngb_mean_abs_u_V_j += (mi / rho_i) * wj * fabsf(rho_i * u_i_prev);
}

/**
 * @brief Band-specific pairwise contribution to particle i's `grad(u)`
 * accumulator (and mirrored contribution to particle j's), the
 * difference-on-(rho*u) form.
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
 * @param grad_u_i (return, accumulated) Particle i's grad(u) accumulator.
 * @param grad_u_j (return, accumulated) Particle j's grad(u) accumulator.
 */
__attribute__((always_inline)) INLINE static void
radiation_gradient_accumulate_band(const float dx[3], float r_inv, float wi_dr,
                                   float wj_dr, float mi, float mj, float rho_i,
                                   float rho_j, float u_i, float u_j,
                                   float grad_u_i[3], float grad_u_j[3]) {

  const float d_ij = rho_i * u_i - rho_j * u_j;
  const float fac_i = -mj * d_ij * wi_dr * r_inv / (rho_i * rho_i);
  const float fac_j = -mi * d_ij * wj_dr * r_inv / (rho_j * rho_j);

  grad_u_i[0] += fac_i * dx[0];
  grad_u_i[1] += fac_i * dx[1];
  grad_u_i[2] += fac_i * dx[2];
  grad_u_j[0] += fac_j * dx[0];
  grad_u_j[1] += fac_j * dx[1];
  grad_u_j[2] += fac_j * dx[2];
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
 * @brief Band-specific pairwise contribution to the Stage-4 anticipatory
 * noise indicator's two kernel-weighted sums over the neighbours' `div(F)`
 * (design-lw-fuv-design-b-dissipation.md Section 5.2), Rosswog 2015a
 * Eq. 87-89 transcribed from `div v` to `div F`.
 *
 * Accumulated in the gradient loop, where `div(F)` is already final. The
 * common `sum_j W_ij` normalisation cancels out of the indicator's own
 * ratio (radiation_isrf.c's #radiation_dissipation_noise_alpha_band), so it
 * is not accumulated.
 *
 * @param wi Particle i's own kernel value, W(r/h_i)*h_i^-dim.
 * @param wj Particle j's own kernel value, W(r/h_j)*h_j^-dim.
 * @param psi_i Particle i's finalized `div(F)` (this band).
 * @param psi_j Particle j's finalized `div(F)` (this band).
 * @param ngb_sum_i (return, accumulated) Particle i's `sum_j W_ij div_F_j`.
 * @param ngb_sum_j (return, accumulated) Particle j's `sum_i W_ji div_F_i`.
 * @param ngb_sum_abs_i (return, accumulated) Particle i's
 * `sum_j W_ij |div_F_j|`.
 * @param ngb_sum_abs_j (return, accumulated) Particle j's
 * `sum_i W_ji |div_F_i|`.
 */
__attribute__((always_inline)) INLINE static void
radiation_noise_indicator_accumulate_band(float wi, float wj, float psi_i,
                                          float psi_j, float *ngb_sum_i,
                                          float *ngb_sum_j,
                                          float *ngb_sum_abs_i,
                                          float *ngb_sum_abs_j) {

  *ngb_sum_i += wi * psi_j;
  *ngb_sum_j += wj * psi_i;
  *ngb_sum_abs_i += wi * fabsf(psi_j);
  *ngb_sum_abs_j += wj * fabsf(psi_i);
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

  radiation_dissipation_accumulate_band(
      dx, r, hi, hj, wi, wj, wi_dr, wj_dr, mi, mj, rho_i, rho_j, fdi->c_hyp,
      fdj->c_hyp, fdi->dissipation_alpha_FUV, fdj->dissipation_alpha_FUV,
      fdi->u_FUV_prev, fdj->u_FUV_prev, g_FUV_i, g_FUV_j,
      &fdi->dissipation_u_FUV, &fdj->dissipation_u_FUV,
      &fdi->ngb_mean_abs_u_V_FUV, &fdj->ngb_mean_abs_u_V_FUV);
  radiation_dissipation_accumulate_band(
      dx, r, hi, hj, wi, wj, wi_dr, wj_dr, mi, mj, rho_i, rho_j, fdi->c_hyp,
      fdj->c_hyp, fdi->dissipation_alpha_LW, fdj->dissipation_alpha_LW,
      fdi->u_LW_prev, fdj->u_LW_prev, g_LW_i, g_LW_j, &fdi->dissipation_u_LW,
      &fdj->dissipation_u_LW, &fdi->ngb_mean_abs_u_V_LW,
      &fdj->ngb_mean_abs_u_V_LW);
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
  float unused_dissipation_u_FUV = 0.f;
  float unused_dissipation_u_LW = 0.f;
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

  radiation_dissipation_accumulate_band(
      dx, r, hi, hj, wi, wj, wi_dr, wj_dr, mi, mj, rho_i, rho_j, fdi->c_hyp,
      fdj->c_hyp, fdi->dissipation_alpha_FUV, fdj->dissipation_alpha_FUV,
      fdi->u_FUV_prev, fdj->u_FUV_prev, g_FUV_i, g_FUV_j,
      &fdi->dissipation_u_FUV, &unused_dissipation_u_FUV,
      &fdi->ngb_mean_abs_u_V_FUV, &unused_ngb_mean_abs_u_V_FUV);
  radiation_dissipation_accumulate_band(
      dx, r, hi, hj, wi, wj, wi_dr, wj_dr, mi, mj, rho_i, rho_j, fdi->c_hyp,
      fdj->c_hyp, fdi->dissipation_alpha_LW, fdj->dissipation_alpha_LW,
      fdi->u_LW_prev, fdj->u_LW_prev, g_LW_i, g_LW_j, &fdi->dissipation_u_LW,
      &unused_dissipation_u_LW, &fdi->ngb_mean_abs_u_V_LW,
      &unused_ngb_mean_abs_u_V_LW);
}

/**
 * @brief `grad(u)` interaction between two particles (symmetric): both
 * particles' accumulators are updated.
 *
 * Runs in the gradient loop, after the density ghost has finalized `u_FUV`/
 * `u_LW` for this step (the exact-relaxation `u` update, radiation_isrf.c):
 * reads them directly, not a `_prev` snapshot, since neither field is
 * written again until star feedback injection, which runs after this loop.
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
  wi *= pow_dimension(hi_inv);
  wj *= pow_dimension(hj_inv);

  struct feedback_part_data *fdi = &pi->feedback_data;
  struct feedback_part_data *fdj = &pj->feedback_data;
  const float rho_i = fdi->rho_prev;
  const float rho_j = fdj->rho_prev;
  const float mi = hydro_get_mass(pi);
  const float mj = hydro_get_mass(pj);

  radiation_gradient_accumulate_band(dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i,
                                     rho_j, fdi->u_FUV, fdj->u_FUV,
                                     fdi->grad_u_FUV, fdj->grad_u_FUV);
  radiation_gradient_accumulate_band(dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i,
                                     rho_j, fdi->u_LW, fdj->u_LW,
                                     fdi->grad_u_LW, fdj->grad_u_LW);

  /* Stages 3 and 4 own per-particle fields that only exist when the stage
     is built, so their call sites are guarded rather than gated on a
     runtime flag; the formulas above stay compiled in either state. */
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

#ifdef RADIATION_LW_FUV_DISSIPATION_ANTICIPATORY_TRIGGER
  radiation_noise_indicator_accumulate_band(
      wi, wj, fdi->div_specific_flux_FUV, fdj->div_specific_flux_FUV,
      &fdi->ngb_sum_div_specific_flux_FUV, &fdj->ngb_sum_div_specific_flux_FUV,
      &fdi->ngb_sum_abs_div_specific_flux_FUV,
      &fdj->ngb_sum_abs_div_specific_flux_FUV);
  radiation_noise_indicator_accumulate_band(
      wi, wj, fdi->div_specific_flux_LW, fdj->div_specific_flux_LW,
      &fdi->ngb_sum_div_specific_flux_LW, &fdj->ngb_sum_div_specific_flux_LW,
      &fdi->ngb_sum_abs_div_specific_flux_LW,
      &fdj->ngb_sum_abs_div_specific_flux_LW);
#else
  (void)wi;
  (void)wj;
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
  wi *= pow_dimension(hi_inv);
  wj *= pow_dimension(hj_inv);

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

  radiation_gradient_accumulate_band(dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i,
                                     rho_j, fdi->u_FUV, fdj->u_FUV,
                                     fdi->grad_u_FUV, unused_grad_u_FUV);
  radiation_gradient_accumulate_band(dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i,
                                     rho_j, fdi->u_LW, fdj->u_LW,
                                     fdi->grad_u_LW, unused_grad_u_LW);

  /* See the symmetric variant above for why these two are guarded rather
     than gated on a runtime flag. */
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

#ifdef RADIATION_LW_FUV_DISSIPATION_ANTICIPATORY_TRIGGER
  float unused_ngb_sum_FUV = 0.f;
  float unused_ngb_sum_LW = 0.f;
  float unused_ngb_sum_abs_FUV = 0.f;
  float unused_ngb_sum_abs_LW = 0.f;

  radiation_noise_indicator_accumulate_band(
      wi, wj, fdi->div_specific_flux_FUV, fdj->div_specific_flux_FUV,
      &fdi->ngb_sum_div_specific_flux_FUV, &unused_ngb_sum_FUV,
      &fdi->ngb_sum_abs_div_specific_flux_FUV, &unused_ngb_sum_abs_FUV);
  radiation_noise_indicator_accumulate_band(
      wi, wj, fdi->div_specific_flux_LW, fdj->div_specific_flux_LW,
      &fdi->ngb_sum_div_specific_flux_LW, &unused_ngb_sum_LW,
      &fdi->ngb_sum_abs_div_specific_flux_LW, &unused_ngb_sum_abs_LW);
#else
  (void)wi;
  (void)wj;
#endif
}

#endif /* SWIFT_RADIATION_PROPAGATION_IACT_GEAR_H */
