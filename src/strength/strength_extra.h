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
#ifndef SWIFT_STRENGTH_EXTRA_H
#define SWIFT_STRENGTH_EXTRA_H

/**
 * @file strength/strength_extra.h
 * @brief Extra strength calculations that are not model-specific, but are only
 * required with certain models.
 */

#include "math.h"
#include "symmetric_matrix.h"

/**
 * @brief Sets the values of additional particle extra properties at a
 * kick time
 *
 * @param p The particle of interest.
 * @param xp The extended data of this particle.
 */
__attribute__((always_inline)) INLINE static void strength_reset_predicted_values_extra(
    struct part *restrict p, const struct xpart *restrict xp) {

#ifdef STRENGTH_DAMAGE_SHEAR_COLLINS
  p->strength_data.strain_tensor = xp->strength_data.strain_tensor_full;
#endif /* STRENGTH_DAMAGE_SHEAR_COLLINS */
}

/**
 * @brief Predict additional particle strength extra properties forward in time when
 * drifting. At beginning of hydro function, before hydro quantities have been drifted.
 *
 * @param p The particle to act upon
 * @param dt_therm The time-step used to evolve hydrodynamical quantities.
 */
__attribute__((always_inline)) INLINE static void strength_predict_extra_beginning(
    struct part *restrict p, const float dt_therm) {

#ifdef STRENGTH_DAMAGE_SHEAR_COLLINS
  float strain_rate_tensor[3][3], rotation_rate_tensor[3][3], rotation_term[3][3];

  /* Compute strain rate and rotation rate. */
  strength_compute_strain_rate_tensor(strain_rate_tensor, p->strength_data.dv_force_loop);
  strength_compute_rotation_rate_tensor(rotation_rate_tensor, p->strength_data.dv_force_loop);

  /* Convert strain tensor to 3x3 float to compute the rotation term . */
  float strain_tensor[3][3];
  get_matrix_from_sym_matrix(strain_tensor,
                             &p->strength_data.strain_tensor);

  /* Compute rotation term. */
  strength_compute_rotation_term(rotation_term, rotation_rate_tensor,
                                 strain_tensor);

  /* Evolve strain tensor. */
  strength_evolve_strain_tensor(strain_tensor, strain_rate_tensor, rotation_term, dt_therm);
  get_sym_matrix_from_matrix(&p->strength_data.strain_tensor, strain_tensor);
#endif /* STRENGTH_DAMAGE_SHEAR_COLLINS */
}

/**
 * @brief Kick the additional particle strength extra properties.
 * At beginning of hydro function, before hydro quantities have been kicked.
 *
 * Additional hydrodynamic quantites are kicked forward in time here. These
 * include thermal quantities (thermal energy or total energy or entropy, ...).
 *
 * @param p The particle to act upon.
 * @param xp The particle extended data to act upon.
 * @param dt_therm The time-step for this kick (for thermodynamic quantities).
 */
__attribute__((always_inline)) INLINE static void strength_kick_extra_beginning(
    struct part *restrict p, struct xpart *restrict xp, const float dt_therm) {

#ifdef STRENGTH_DAMAGE_SHEAR_COLLINS
  float strain_rate_tensor[3][3], rotation_rate_tensor[3][3], rotation_term[3][3];

  /* Compute strain rate and rotation rate. */
  strength_compute_strain_rate_tensor(strain_rate_tensor, p->strength_data.dv_force_loop);
  strength_compute_rotation_rate_tensor(rotation_rate_tensor, p->strength_data.dv_force_loop);

  /* Convert strain tensor to 3x3 float to compute the rotation term . */
  float strain_tensor[3][3];
  get_matrix_from_sym_matrix(strain_tensor,
                             &xp->strength_data.strain_tensor_full);

  /* Compute rotation term. */
  strength_compute_rotation_term(rotation_term, rotation_rate_tensor,
                                 strain_tensor);

  /* Evolve strain tensor. */
  strength_evolve_strain_tensor(strain_tensor, strain_rate_tensor, rotation_term, dt_therm);
  get_sym_matrix_from_matrix(&xp->strength_data.strain_tensor_full, strain_tensor);
#endif /* STRENGTH_DAMAGE_SHEAR_COLLINS */
}

/**
 * @brief Initialises the extra stress properties for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions or assignments between the particle
 * and extended particle fields.
 *
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static void strength_first_init_part_extra(
    struct part *restrict p, struct xpart *restrict xp) {

#ifdef STRENGTH_DAMAGE_SHEAR_COLLINS
  zero_sym_matrix(&p->strength_data.strain_tensor);
  zero_sym_matrix(&xp->strength_data.strain_tensor_full);
#endif /* STRENGTH_DAMAGE_SHEAR_COLLINS */
}

#endif /* SWIFT_STRENGTH_EXTRA_H */
