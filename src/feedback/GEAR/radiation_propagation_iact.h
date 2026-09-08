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
  (void)wi;
  (void)wj;
  const float wi_dr = wi_dx * pow_dimension_plus_one(hi_inv);
  const float wj_dr = wj_dx * pow_dimension_plus_one(hj_inv);

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
  (void)wi;
  (void)wj;
  const float wi_dr = wi_dx * pow_dimension_plus_one(hi_inv);
  const float wj_dr = wj_dx * pow_dimension_plus_one(hj_inv);

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

  radiation_divergence_accumulate_band(
      dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, fdi->specific_flux_FUV,
      fdj->specific_flux_FUV, &fdi->div_specific_flux_FUV,
      &unused_div_specific_flux_FUV);
  radiation_divergence_accumulate_band(
      dx, r_inv, wi_dr, wj_dr, mi, mj, rho_i, rho_j, fdi->specific_flux_LW,
      fdj->specific_flux_LW, &fdi->div_specific_flux_LW,
      &unused_div_specific_flux_LW);
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
  (void)wi;
  (void)wj;
  const float wi_dr = wi_dx * pow_dimension_plus_one(hi_inv);
  const float wj_dr = wj_dx * pow_dimension_plus_one(hj_inv);

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
  (void)wi;
  (void)wj;
  const float wi_dr = wi_dx * pow_dimension_plus_one(hi_inv);
  const float wj_dr = wj_dx * pow_dimension_plus_one(hj_inv);

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
}

#endif /* SWIFT_RADIATION_PROPAGATION_IACT_GEAR_H */
