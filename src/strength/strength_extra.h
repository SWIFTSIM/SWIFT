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
  p->strength_data.total_plastic_strain = xp->strength_data.total_plastic_strain_full;
#endif /* STRENGTH_DAMAGE_SHEAR_COLLINS */
}

/**
 * @brief Predict additional particle strength extra properties forward in time when
 * drifting. At beginning of hydro function, before hydro quantities have been drifted.
 *
 * @param p The particle to act upon
 * @param density The density.
 * @param u The specific internal energy.
 * @param yield_stress The yield stress.
 * @param dt_therm The time-step used to evolve hydrodynamical quantities.
 */
__attribute__((always_inline)) INLINE static void strength_predict_extra_beginning(
    struct part *restrict p, const float density, const float u, const float yield_stress,  const float dt_therm) {

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

  /* Make copy of strain tensor prior to apply yield criterion. */
  const struct sym_matrix strain_tensor_pre_yield = p->strength_data.strain_tensor;

  /* Apply yield stress to strain tensor. */
  yield_apply_yield_stress_to_sym_matrix(
         &p->strength_data.strain_tensor, p->strength_data.deviatoric_stress_tensor, density, u, yield_stress);

  /* Calculate the change in elastic strain in the time-step. */
  struct sym_matrix delta_elastic_strain_tensor;
  delta_elastic_strain_tensor.xx = p->strength_data.strain_tensor.xx - strain_tensor_pre_yield.xx;
  delta_elastic_strain_tensor.yy = p->strength_data.strain_tensor.yy - strain_tensor_pre_yield.yy;
  delta_elastic_strain_tensor.zz = p->strength_data.strain_tensor.zz - strain_tensor_pre_yield.zz;
  delta_elastic_strain_tensor.xy = p->strength_data.strain_tensor.xy - strain_tensor_pre_yield.xy;
  delta_elastic_strain_tensor.xz = p->strength_data.strain_tensor.xz - strain_tensor_pre_yield.xz;
  delta_elastic_strain_tensor.yz = p->strength_data.strain_tensor.yz - strain_tensor_pre_yield.yz;

  /* Add contribution to measure of total plastic strain. */
  p->strength_data.total_plastic_strain += strength_compute_sym_matrix_J_2(delta_elastic_strain_tensor);
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
 * @param density The density.
 * @param u The specific internal energy.
 * @param yield_stress The yield stress.
 * @param dt_therm The time-step for this kick (for thermodynamic quantities).
 */
__attribute__((always_inline)) INLINE static void strength_kick_extra_beginning(
    struct part *restrict p, struct xpart *restrict xp,  const float density, const float u, const float yield_stress, const float dt_therm) {

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

  /* Make copy of strain tensor prior to apply yield criterion. */
  const struct sym_matrix strain_tensor_pre_yield = xp->strength_data.strain_tensor_full;

  /* Apply yield stress to strain tensor. */
  yield_apply_yield_stress_to_sym_matrix(
         &xp->strength_data.strain_tensor_full, xp->strength_data.deviatoric_stress_tensor_full, density, u, yield_stress);

 /* Calculate the change in elastic strain in the time-step. */
 struct sym_matrix delta_elastic_strain_tensor;
 delta_elastic_strain_tensor.xx = xp->strength_data.strain_tensor_full.xx - strain_tensor_pre_yield.xx;
 delta_elastic_strain_tensor.yy = xp->strength_data.strain_tensor_full.yy - strain_tensor_pre_yield.yy;
 delta_elastic_strain_tensor.zz = xp->strength_data.strain_tensor_full.zz - strain_tensor_pre_yield.zz;
 delta_elastic_strain_tensor.xy = xp->strength_data.strain_tensor_full.xy - strain_tensor_pre_yield.xy;
 delta_elastic_strain_tensor.xz = xp->strength_data.strain_tensor_full.xz - strain_tensor_pre_yield.xz;
 delta_elastic_strain_tensor.yz = xp->strength_data.strain_tensor_full.yz - strain_tensor_pre_yield.yz;

 /* Add contribution to measure of total plastic strain. */
 xp->strength_data.total_plastic_strain_full += strength_compute_sym_matrix_J_2(delta_elastic_strain_tensor);
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

  p->strength_data.total_plastic_strain = 0.f;
  xp->strength_data.total_plastic_strain_full = 0.f;
#endif /* STRENGTH_DAMAGE_SHEAR_COLLINS */
}

#endif /* SWIFT_STRENGTH_EXTRA_H */
