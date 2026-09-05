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
 * e^{-tau}-weighted kernel sums needed for the normalized mixing weight;
 * feedback_end_density() consumes them once per particle after the
 * density loop completes.
 */

#include "chemistry.h"
#include "hydro.h"
#include "kernel_hydro.h"
#include "radiation.h"

/**
 * @brief Band-specific pairwise contribution to particle i's propagation
 * accumulators from particle j's current field value.
 *
 * @param r Comoving particle separation.
 * @param wi Kernel weight evaluated at particle i's own smoothing length.
 * @param Z_i Particle i's metal mass fraction.
 * @param rho_i_phys Particle i's physical density.
 * @param Z_j Particle j's metal mass fraction.
 * @param rho_j_phys Particle j's physical density.
 * @param u_j Particle j's current field value (u_FUV or u_LW).
 * @param sigma_d_band_cgs Band cross-section per hydrogen nucleon, cm^2.
 * @param us Unit system.
 * @param sum_w (return, accumulated) Weighted normalization denominator.
 * @param sum_wu (return, accumulated) Weighted field numerator.
 */
__attribute__((always_inline)) INLINE static void
radiation_propagation_accumulate_band(
    float r, float wi, float Z_i, float rho_i_phys, float Z_j,
    float rho_j_phys, float u_j, float sigma_d_band_cgs,
    const struct unit_system *us, float *sum_w, float *sum_wu) {

  const float kappa_i =
      radiation_get_part_linear_absorption_rate(us, Z_i, rho_i_phys,
                                                sigma_d_band_cgs);
  const float kappa_j =
      radiation_get_part_linear_absorption_rate(us, Z_j, rho_j_phys,
                                                sigma_d_band_cgs);
  const float kappa_sum = kappa_i + kappa_j;
  const float kappa_ij =
      kappa_sum > 0.0f ? 2.0f * kappa_i * kappa_j / kappa_sum : 0.0f;
  const float tau_ij = kappa_ij * r;
  const float w = wi * expf(-tau_ij);

  *sum_w += w;
  *sum_wu += w * u_j;
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
 * @param us Unit system.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_isrf_propagation(const float r2, const float dx[3],
                             const float hi, const float hj,
                             struct part *restrict pi,
                             struct part *restrict pj, const float a,
                             const float H, const struct unit_system *us) {

  const float r = sqrtf(r2);
  const float a3_inv = 1.0f / (a * a * a);
  const float Z_i = chemistry_get_total_metal_mass_fraction_for_cooling(pi);
  const float Z_j = chemistry_get_total_metal_mass_fraction_for_cooling(pj);
  const float rho_i_phys = pi->rho * a3_inv;
  const float rho_j_phys = pj->rho * a3_inv;

  float wi, wj;
  kernel_eval(r / hi, &wi);
  kernel_eval(r / hj, &wj);

  radiation_propagation_accumulate_band(
      r, wi, Z_i, rho_i_phys, Z_j, rho_j_phys, pj->feedback_data.u_FUV,
      RADIATION_SIGMA_D_FUV_CGS, us, &pi->feedback_data.isrf_prop_sum_w_FUV,
      &pi->feedback_data.isrf_prop_sum_wu_FUV);
  radiation_propagation_accumulate_band(
      r, wi, Z_i, rho_i_phys, Z_j, rho_j_phys, pj->feedback_data.u_LW,
      RADIATION_SIGMA_D_LW_CGS, us, &pi->feedback_data.isrf_prop_sum_w_LW,
      &pi->feedback_data.isrf_prop_sum_wu_LW);

  radiation_propagation_accumulate_band(
      r, wj, Z_j, rho_j_phys, Z_i, rho_i_phys, pi->feedback_data.u_FUV,
      RADIATION_SIGMA_D_FUV_CGS, us, &pj->feedback_data.isrf_prop_sum_w_FUV,
      &pj->feedback_data.isrf_prop_sum_wu_FUV);
  radiation_propagation_accumulate_band(
      r, wj, Z_j, rho_j_phys, Z_i, rho_i_phys, pi->feedback_data.u_LW,
      RADIATION_SIGMA_D_LW_CGS, us, &pj->feedback_data.isrf_prop_sum_w_LW,
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
 * @param us Unit system.
 */
__attribute__((always_inline)) INLINE static void
runner_iact_nonsym_isrf_propagation(const float r2, const float dx[3],
                                    const float hi, const float hj,
                                    struct part *restrict pi,
                                    const struct part *restrict pj,
                                    const float a, const float H,
                                    const struct unit_system *us) {

  const float r = sqrtf(r2);
  const float a3_inv = 1.0f / (a * a * a);
  const float Z_i = chemistry_get_total_metal_mass_fraction_for_cooling(pi);
  const float Z_j = chemistry_get_total_metal_mass_fraction_for_cooling(pj);
  const float rho_i_phys = pi->rho * a3_inv;
  const float rho_j_phys = pj->rho * a3_inv;

  float wi;
  kernel_eval(r / hi, &wi);

  radiation_propagation_accumulate_band(
      r, wi, Z_i, rho_i_phys, Z_j, rho_j_phys, pj->feedback_data.u_FUV,
      RADIATION_SIGMA_D_FUV_CGS, us, &pi->feedback_data.isrf_prop_sum_w_FUV,
      &pi->feedback_data.isrf_prop_sum_wu_FUV);
  radiation_propagation_accumulate_band(
      r, wi, Z_i, rho_i_phys, Z_j, rho_j_phys, pj->feedback_data.u_LW,
      RADIATION_SIGMA_D_LW_CGS, us, &pi->feedback_data.isrf_prop_sum_w_LW,
      &pi->feedback_data.isrf_prop_sum_wu_LW);
}

#endif /* SWIFT_RADIATION_PROPAGATION_IACT_GEAR_H */
