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
 * @brief Gas-gas density-loop hook for the Yukawa screened-diffusion
 * propagation of the u_FUV/u_LW fields, mirroring GEAR chemistry's own
 * smoothed-metallicity diffusion (chemistry_iact.h). Accumulates the
 * e^{-tau}-weighted kernel sums needed for the normalized mixing weight,
 * reading each neighbour's stable u_*_prev snapshot (not the live,
 * per-h-iteration u_FUV/u_LW output) so the result does not depend on
 * h-iteration count or convergence order; feedback_end_density() consumes
 * the sums once per h-iteration.
 */

#include "kernel_hydro.h"
#include "radiation.h"

/**
 * @brief Band-specific pairwise contribution to particle i's propagation
 * accumulators from particle j's stable field snapshot.
 *
 * @param r Comoving particle separation.
 * @param wi Kernel weight evaluated at particle i's own smoothing length.
 * @param kappa_i Particle i's cached absorption rate for this band.
 * @param kappa_j Particle j's cached absorption rate for this band.
 * @param u_j_prev Particle j's stable field snapshot (u_FUV_prev/u_LW_prev).
 * @param sum_w (return, accumulated) Weighted normalization denominator.
 * @param sum_wu (return, accumulated) Weighted field numerator.
 */
__attribute__((always_inline)) INLINE static void
radiation_propagation_accumulate_band(float r, float wi, float kappa_i,
                                      float kappa_j, float u_j_prev,
                                      float *sum_w, float *sum_wu) {

  const float kappa_sum = kappa_i + kappa_j;
  const float kappa_ij =
      kappa_sum > 0.0f ? 2.0f * kappa_i * kappa_j / kappa_sum : 0.0f;
  const float w = wi * expf(-kappa_ij * r);

  *sum_w += w;
  *sum_wu += w * u_j_prev;
}

/**
 * @brief Yukawa propagation interaction between two particles (symmetric).
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle.
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param us Unit system (unused: absorption rates are precomputed and
 * cached by radiation_init_part_propagation).
 */
__attribute__((always_inline)) INLINE static void runner_iact_isrf_propagation(
    const float r2, const float dx[3], const float hi, const float hj,
    struct part *restrict pi, struct part *restrict pj, const float a,
    const float H, const struct unit_system *us) {

  const float r = sqrtf(r2);
  float wi, wj;
  kernel_eval(r / hi, &wi);
  kernel_eval(r / hj, &wj);

  radiation_propagation_accumulate_band(
      r, wi, pi->feedback_data.kappa_FUV, pj->feedback_data.kappa_FUV,
      pj->feedback_data.u_FUV_prev, &pi->feedback_data.isrf_prop_sum_w_FUV,
      &pi->feedback_data.isrf_prop_sum_wu_FUV);
  radiation_propagation_accumulate_band(
      r, wi, pi->feedback_data.kappa_LW, pj->feedback_data.kappa_LW,
      pj->feedback_data.u_LW_prev, &pi->feedback_data.isrf_prop_sum_w_LW,
      &pi->feedback_data.isrf_prop_sum_wu_LW);

  radiation_propagation_accumulate_band(
      r, wj, pj->feedback_data.kappa_FUV, pi->feedback_data.kappa_FUV,
      pi->feedback_data.u_FUV_prev, &pj->feedback_data.isrf_prop_sum_w_FUV,
      &pj->feedback_data.isrf_prop_sum_wu_FUV);
  radiation_propagation_accumulate_band(
      r, wj, pj->feedback_data.kappa_LW, pi->feedback_data.kappa_LW,
      pi->feedback_data.u_LW_prev, &pj->feedback_data.isrf_prop_sum_w_LW,
      &pj->feedback_data.isrf_prop_sum_wu_LW);
}

/**
 * @brief Yukawa propagation interaction between two particles
 * (non-symmetric): only particle i's accumulators are updated.
 *
 * @param r2 Comoving square distance between the two particles.
 * @param dx Comoving vector separating both particles (pi - pj).
 * @param hi Comoving smoothing-length of particle i.
 * @param hj Comoving smoothing-length of particle j.
 * @param pi First particle.
 * @param pj Second particle (not updated).
 * @param a Current scale factor.
 * @param H Current Hubble parameter.
 * @param us Unit system (unused, see runner_iact_isrf_propagation).
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_isrf_propagation(const float r2, const float dx[3],
                                    const float hi, const float hj,
                                    struct part *restrict pi,
                                    const struct part *restrict pj,
                                    const float a, const float H,
                                    const struct unit_system *us) {

  const float r = sqrtf(r2);
  float wi;
  kernel_eval(r / hi, &wi);

  radiation_propagation_accumulate_band(
      r, wi, pi->feedback_data.kappa_FUV, pj->feedback_data.kappa_FUV,
      pj->feedback_data.u_FUV_prev, &pi->feedback_data.isrf_prop_sum_w_FUV,
      &pi->feedback_data.isrf_prop_sum_wu_FUV);
  radiation_propagation_accumulate_band(
      r, wi, pi->feedback_data.kappa_LW, pj->feedback_data.kappa_LW,
      pj->feedback_data.u_LW_prev, &pi->feedback_data.isrf_prop_sum_w_LW,
      &pi->feedback_data.isrf_prop_sum_wu_LW);
}

#endif /* SWIFT_RADIATION_PROPAGATION_IACT_GEAR_H */
