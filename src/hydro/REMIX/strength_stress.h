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
#include "hydro_kernels.h"
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

  for (int i = 0; i < 6; i++) {
    p->stress_tensor.elements[i] = p->deviatoric_stress_tensor.elements[i];
  }

  float effective_pressure =
      effective_pressure_from_damage(p, pressure);

  p->stress_tensor.xx -= effective_pressure;
  p->stress_tensor.yy -= effective_pressure;
  p->stress_tensor.zz -= effective_pressure;

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
      const float delta_p = 0.5f * (powf(pi->mass /
     pi->rho, 1.f/hydro_dimension) + powf(pj->mass /
     pj->rho, 1.f/hydro_dimension));

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
     artif_stress_epsilon * pairwise_stress_tensor_i[i][j]; if
     (pairwise_stress_tensor_j[i][j] > 0.f) pairwise_stress_tensor_j[i][j] -=
     artif_stress_f * artif_stress_epsilon * pairwise_stress_tensor_j[i][j];
          }
      }
   #elif defined(STRENGTH_STRESS_BASIS_INDP)
  /* New version that is independent of basis:
   * Principal stresses are basis-independent
   * Adding the same constant to all diagonal elemenets is equivalent in any
   * basis: M + c*I = P (M + c*I) P^-1 = P M P^-1 + P (c*I) P^-1 = P M P^-1 +
   * c*I This is a lot easier than e.g. Gray, Monaghan, and Swift (2001) */
  const float mean_h = 0.5f * (pi->h + pj->h);
  const float delta_p =
      0.5f * (powf(pi->mass / pi->rho, 1.f / hydro_dimension) +
              powf(pj->mass / pj->rho, 1.f / hydro_dimension));

  float wij_delta_p;
  kernel_eval(delta_p / mean_h, &wij_delta_p);

  float wij_r;
  kernel_eval(r / mean_h, &wij_r);

  // This factor should be set in extra parameter file
  const float artif_stress_n = method_artif_stress_n();
  const float artif_stress_f = powf(wij_r / wij_delta_p, artif_stress_n);

  // This factor should be set in extra parameter file
  const float artif_stress_epsilon = method_artif_stress_epsilon();

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
          artif_stress_f * artif_stress_epsilon * max_principal_stress_i;
    if (max_principal_stress_j > 0)
      pairwise_stress_tensor_j[i][i] -=
          artif_stress_f * artif_stress_epsilon * max_principal_stress_j;
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
 * @brief Evolve the deviatoric stress tensor
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void evolve_deviatoric_stress(
    struct part *restrict p, const float dt_therm, const float density,
    const float pressure) {

  if (p->phase_state == mat_phase_state_fluid) {
    // No stress for fluids
    zero_sym_matrix(&p->deviatoric_stress_tensor);
  } else {
    // Update solid stress
    for (int i = 0; i < 6; i++) {
      p->deviatoric_stress_tensor.elements[i] +=
          p->dS_dt.elements[i] * dt_therm;
    }
    // Yield stress
    adjust_deviatoric_stress_tensor_by_yield_stress(
        p, &p->deviatoric_stress_tensor, density, pressure);
  }
}

#endif /* MATERIAL_STRENGTH */
#endif /* SWIFT_PLANETARY_STRESS_H */
