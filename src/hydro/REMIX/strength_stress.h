/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
#ifndef SWIFT_PLANETARY_STRESS_H
#define SWIFT_PLANETARY_STRESS_H
#ifdef MATERIAL_STRENGTH

/**
 * @file Planetary/strength_stress.h
 * @brief REMIX implementation of SPH with material strength
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"
#include "strength_damage.h"
#include "strength_utilities.h"
#include "strength_yield.h"

/**
 * @brief Set the (symmetric) stress tensor by combining the deviatoric with the
 * pressure.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_set_stress_tensor(
    struct part *restrict p, const float pressure) {

  #ifdef STRENGTH_DAMAGE
    p->stress_tensor = stress_tensor_damaged(p->deviatoric_stress_tensor, pressure, p->damage);
  #else
    p->stress_tensor = p->deviatoric_stress_tensor;
    p->stress_tensor.xx -= pressure;
    p->stress_tensor.yy -= pressure;
    p->stress_tensor.zz -= pressure;
  #endif /* STRENGTH_DAMAGE */  
    
    // Compute principal stresses
    sym_matrix_compute_eigenvalues(p->principal_stress_eigen, p->stress_tensor);
}

/**
 * @brief Update the pairwise stress tensors with artificial stress.
 */
__attribute__((always_inline)) INLINE static void strength_add_artif_stress(
    float pairwise_stress_tensor_i[3][3], float pairwise_stress_tensor_j[3][3],
    const struct part *restrict pi, const struct part *restrict pj,
    const float r) {
   #if defined(STRENGTH_STRESS_MON2000)
      // Artificial stress (Monaghan, 2000)
      const float mean_h = 0.5f * (pi->h + pj->h);
      const float delta_p = mean_h / 1.487; // ### hardcoded for now

      float wij_delta_p;
      kernel_eval(delta_p / mean_h, &wij_delta_p);

      float wij_r;
      kernel_eval(r / mean_h, &wij_r);

      // This factor should be set in extra parameter file
      const float artif_stress_n = method_artif_stress_n();
      const float artif_stress_f = powf(wij_r / wij_delta_p, artif_stress_n);

      // This factor should be set in extra parameter file
      const float artif_stress_epsilon = method_artif_stress_epsilon();

      for (int i = 0; i < 3; ++i) {
          for (int j = 0; j < 3; ++j) {
            if (pairwise_stress_tensor_i[i][j] > 0.f)
              pairwise_stress_tensor_i[i][j] -= artif_stress_f *
                  artif_stress_epsilon * pairwise_stress_tensor_i[i][j]; 
            if (pairwise_stress_tensor_j[i][j] > 0.f) 
              pairwise_stress_tensor_j[i][j] -= artif_stress_f * 
                  artif_stress_epsilon * pairwise_stress_tensor_j[i][j];
          }
      }
   #elif defined(STRENGTH_STRESS_BASIS_INDP)
  /* New version that is independent of basis:
   * Principal stresses are basis-independent
   * Adding the same constant to all diagonal elemenets is equivalent in any
   * basis: M + c*I = P (M + c*I) P^-1 = P M P^-1 + P (c*I) P^-1 = P M P^-1 +
   * c*I This is a lot easier than e.g. Gray, Monaghan, and Swift (2001) */

  // ### hardcoded for now so it works with Planetary (WC2)
  const float eta_crit = 1.f / 1.487;//viscosity_global.eta_crit;
  float eta_ab = min(r / pi->h, r / pj->h);

  float artif_stress_f = 0.f;
  if (eta_ab < eta_crit) {
    // ### hardcoded for now so it works with Planetary    
    artif_stress_f = 1.f - expf(-(eta_ab - eta_crit) * (eta_ab - eta_crit) /
                   0.04f);// hydro_slope_limiter_exp_denom);
  }

  float max_principal_stress_i = pi->principal_stress_eigen[0];
  if (pi->principal_stress_eigen[1] > max_principal_stress_i)
    max_principal_stress_i = pi->principal_stress_eigen[1];
  if (pi->principal_stress_eigen[2] > max_principal_stress_i)
    max_principal_stress_i = pi->principal_stress_eigen[2];

  float max_principal_stress_j = pj->principal_stress_eigen[0];
  if (pj->principal_stress_eigen[1] > max_principal_stress_j)
    max_principal_stress_j = pj->principal_stress_eigen[1];
  if (pj->principal_stress_eigen[2] > max_principal_stress_j)
    max_principal_stress_j = pj->principal_stress_eigen[2];

  for (int i = 0; i < 3; ++i) {
    if (max_principal_stress_i > 0)
      pairwise_stress_tensor_i[i][i] -=
          artif_stress_f * max_principal_stress_i;
    if (max_principal_stress_j > 0)
      pairwise_stress_tensor_j[i][i] -=
          artif_stress_f * max_principal_stress_j;
  }
  #endif
}

/**
 * @brief Calculates the stress tensor with strength for the force interaction.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void
hydro_set_pairwise_stress_tensors_strength(float pairwise_stress_tensor_i[3][3],
                                           float pairwise_stress_tensor_j[3][3],
                                           const struct part *restrict pi,
                                           const struct part *restrict pj,
                                           const float r) {

  // Use the full stress tensor for solid particles
  if ((pi->phase_state == mat_phase_state_solid) &&
      (pj->phase_state == mat_phase_state_solid)) {

    get_matrix_from_sym_matrix(pairwise_stress_tensor_i, &pi->stress_tensor);
    get_matrix_from_sym_matrix(pairwise_stress_tensor_j, &pj->stress_tensor);

    strength_add_artif_stress(pairwise_stress_tensor_i,
                              pairwise_stress_tensor_j, pi, pj, r);
  }
}

/**
 * @brief Calculate time derivative of the deviatoric stress tensor
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void calculate_dS_dt(struct part *restrict p) {

  float strain_rate_tensor[3][3], rotation_rate_tensor[3][3],
      rotation_term[3][3];
  float deviatoric_stress_tensor[3][3];
  get_matrix_from_sym_matrix(deviatoric_stress_tensor,
                             &p->deviatoric_stress_tensor);

  // Set the strain and rotation rates
  calculate_strain_rate_tensor(p, strain_rate_tensor);
  calculate_rotation_rate_tensor(p, rotation_rate_tensor);
  calculate_rotation_term(rotation_term, rotation_rate_tensor,
                          deviatoric_stress_tensor);  
    
  // Compute time derivative of the deviatoric stress tensor (Hooke's law)
  const float shear_mod = material_shear_mod(p->mat_id);
  float dS_dt[3][3];
  for (int i = 0; i < 3; i++) {
    for (int j = 0; j < 3; j++) {
      dS_dt[i][j] =
          2.0f * shear_mod * strain_rate_tensor[i][j] + rotation_term[i][j];
    }

    dS_dt[i][i] -= 2.0f * shear_mod *
                   (strain_rate_tensor[0][0] + strain_rate_tensor[1][1] +
                    strain_rate_tensor[2][2]) /
                   3.f;
  }

  get_sym_matrix_from_matrix(&p->dS_dt, dS_dt);
}

/**
 * @brief Evolve the deviatoric stress tensor
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void evolve_deviatoric_stress(
    struct part *restrict p, struct sym_matrix *deviatoric_stress_tensor, const int phase_state, 
    float dt_therm) {

  if (phase_state == mat_phase_state_fluid) {
    // No stress for fluids
    zero_sym_matrix(deviatoric_stress_tensor);
  } else {
    // Update solid stress
    for (int i = 0; i < 6; i++) {
      deviatoric_stress_tensor->elements[i] +=
          p->dS_dt.elements[i] * dt_therm;
    }
  }    
}

#endif /* MATERIAL_STRENGTH */
#endif /* SWIFT_PLANETARY_STRESS_H */
