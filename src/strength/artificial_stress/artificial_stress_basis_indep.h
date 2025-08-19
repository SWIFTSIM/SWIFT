/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2025 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
 * @brief Basis independent artificial stress method.
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Apply artificial stress to pairwise stress tensors.
 *
 * Adds a negative factor to diagonal elements of the stress tensor if any of
 * the principal stresses are positive. The same factor is used for each of the
 * diagonal elements, and it consists of the maximum principal stress multiplied
 * by a Gaussian function that increases from 0 to 1 with reduced particle
 * separation. Similar separation functions  are included in the artificial
 * viscosity and diffusion slope limiters in REMIX. Adding the same constant to
 * all diagonal elemenets is equivalent in any basis:
 * M + c*I = P (M + c*I) P^-1 = P M P^-1 + P (c*I) P^-1 = P M P^-1 + c*I
 * Unlike the basis independent method presented by Gray+2001, this method does
 * not require basis transformation of the stress tensor, and it is therefore
 * more simple and less computationally intesive.
 *
 * @param pairwise_stress_tensor_i Stress tensor of particle i for its interactiion with j.
 * @param pairwise_stress_tensor_j Stress tensor of particle j for its interactiion with i.
 * @param pi First particle.
 * @param pj Second particle.
 * @param r The particle separation.
 */
__attribute__((always_inline)) INLINE static void artif_stress_apply_artif_stress_to_pairwise_stress_tensors(
    float pairwise_stress_tensor_i[3][3], float pairwise_stress_tensor_j[3][3],
    const struct part *restrict pi, const struct part *restrict pj,
    const float r) {

 /* Calculate qunatities needed for separation factor. */
  const float eta_crit = 1.f / 1.487;//viscosity_global.eta_crit; // ### hardcoded for now so it works with Planetary (WC2)
  const float eta_ab = fminf(r / pi->h, r / pj->h);

  /* No artificial stress is applied for scaled particle separations greater
   * than eta_crit. */
  if (eta_ab >= eta_crit) {
    return;
  }

  /* Calculate separation factor, artif_stress_f. */
  const float artif_stress_f = 1.f - expf(-(eta_ab - eta_crit) * (eta_ab - eta_crit) /
                   0.04f);// ### hardcoded for now hydro_slope_limiter_exp_denom;

  /* Get max principal stresses for particles i and j. */
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

  /* Apply artificial stress. */
  if (max_principal_stress_i > 0.f) {
    pairwise_stress_tensor_i[0][0] -= artif_stress_f * max_principal_stress_i;
    pairwise_stress_tensor_i[1][1] -= artif_stress_f * max_principal_stress_i;
    pairwise_stress_tensor_i[2][2] -= artif_stress_f * max_principal_stress_i;
  }

  if (max_principal_stress_j > 0.f) {
    pairwise_stress_tensor_j[0][0] -= artif_stress_f * max_principal_stress_j;
    pairwise_stress_tensor_j[1][1] -= artif_stress_f * max_principal_stress_j;
    pairwise_stress_tensor_j[2][2] -= artif_stress_f * max_principal_stress_j;
  }
}

#endif /* SWIFT_ARTIFICIAL_STRESS_BASIS_INDEP_H */
