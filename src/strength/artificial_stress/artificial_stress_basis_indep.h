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
#ifndef SWIFT_ARTIFICIAL_STRESS_BASIS_INDEP_H
#define SWIFT_ARTIFICIAL_STRESS_BASIS_INDEP_H

/**
 * @file strength/artificial_stress/artificial_stress_basis_indep.h
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Update the pairwise stress tensors with artificial stress.
 */
__attribute__((always_inline)) INLINE static void strength_add_artif_stress(
    float pairwise_stress_tensor_i[3][3], float pairwise_stress_tensor_j[3][3],
    const struct part *restrict pi, const struct part *restrict pj,
    const float r) {

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

  float max_principal_stress_i = pi->strength_data.principal_stress_eigen[0];
  if (pi->strength_data.principal_stress_eigen[1] > max_principal_stress_i)
    max_principal_stress_i = pi->strength_data.principal_stress_eigen[1];
  if (pi->strength_data.principal_stress_eigen[2] > max_principal_stress_i)
    max_principal_stress_i = pi->strength_data.principal_stress_eigen[2];

  float max_principal_stress_j = pj->strength_data.principal_stress_eigen[0];
  if (pj->strength_data.principal_stress_eigen[1] > max_principal_stress_j)
    max_principal_stress_j = pj->strength_data.principal_stress_eigen[1];
  if (pj->strength_data.principal_stress_eigen[2] > max_principal_stress_j)
    max_principal_stress_j = pj->strength_data.principal_stress_eigen[2];

  for (int i = 0; i < 3; ++i) {
    if (max_principal_stress_i > 0)
      pairwise_stress_tensor_i[i][i] -=
          artif_stress_f * max_principal_stress_i;
    if (max_principal_stress_j > 0)
      pairwise_stress_tensor_j[i][i] -=
          artif_stress_f * max_principal_stress_j;
  }
}

#endif /* SWIFT_ARTIFICIAL_STRESS_BASIS_INDEP_H */