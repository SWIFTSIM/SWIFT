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
#ifndef SWIFT_STRENGTH_STRESS_TENSOR_H
#define SWIFT_STRENGTH_STRESS_TENSOR_H

/**
 * @file strength/strength_stress_tensor.h
 * @brief Hooke's Law model for elastic stress
 */

#include "const.h"
#include "math.h"

/**
 * @brief Updates the max wave speed of solids with the longitudinal wave speed.
 *
 * @param wave_speed The wave speed to be updated.
 * @param p The particle of interest.
 * @param soundspeed The sound speed.
 * @param density The density.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void
strength_compute_max_wave_speed_stress_tensor(float *wave_speed, const struct part *restrict p,
                                              const float soundspeed, const float density,
                                              const struct mat_params *mp) {
  if (p->phase == mat_phase_solid) {
    /* Speed of longitudinal elastic wave. */
    const float shear_mod = mp->shear_mod;
    *wave_speed = sqrtf(soundspeed * soundspeed + (4.f / 3.f) * shear_mod / density);
  }
}

/**
 * @brief Modify the pressure when in tension.
 *
 * Damage scales negative pressures by the damage factor of the yield stress method.
 *
 * @param pressure The pressure, to be modified if negative.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void strength_modify_pressure_in_tension(
    float *pressure, struct strength_state *s, const struct mat_params *mp) {

  if (*pressure < 0.f) {
    *pressure *= yield_damage_factor(s, mp);
  }
}

/**
 * @brief Calculates the pairwise stress tensors for the force interaction.
 *
 * The stress tensors used for the force interaction between a specific pair of
 * particles. These differ from the particle's own stress tensor, since they
 * factor in e.g. the phases of the two particles as well as the contribution of
 * artificial stress for the pairwise interaction.
 *
 * @param pairwise_stress_tensor_i Stress tensor of particle i for its interaction with j.
 * @param pairwise_stress_tensor_j Stress tensor of particle j for its interaction with i.
 * @param pi First particle.
 * @param pj Second particle.
 * @param r The particle separation.
 * @param mp_i The material parameters of the first particle.
 * @param mp_j The material parameters of the second particle.
 * @param strength_props The material strength parameters.
 */
__attribute__((always_inline)) INLINE static void
strength_set_pairwise_stress_tensors(float pairwise_stress_tensor_i[3][3],
                                     float pairwise_stress_tensor_j[3][3],
                                     const struct part *restrict pi,
                                     const struct part *restrict pj,
                                     const float r,
                                     const struct mat_params *mp_i,
                                     const struct mat_params *mp_j,
                                     const struct strength_props *strength_props) {

  /* Only overwrite the fluid pairwise stress tensors for interactions with strength. */
  if (!strength_is_strengthless_interaction(pi, pj)) {

    /* Get stress tensors. */
    get_matrix_from_sym_matrix(pairwise_stress_tensor_i, &pi->strength_data.stress_tensor);
    get_matrix_from_sym_matrix(pairwise_stress_tensor_j, &pj->strength_data.stress_tensor);

    /* Apply artificial stress. */
    artif_stress_apply_artif_stress_to_pairwise_stress_tensors(pairwise_stress_tensor_i,
                                                               pairwise_stress_tensor_j, pi, pj, r,
                                                               mp_i, mp_j, strength_props);
  }
}


/**
 * @brief Sets the values of additional particle stress tensor properties at a
 * kick time
 *
 * @param p The particle of interest.
 * @param xp The extended data of this particle.
 */
__attribute__((always_inline)) INLINE static void strength_reset_predicted_values_stress_tensor(
    struct part *restrict p, const struct xpart *restrict xp) {

  p->strength_data.deviatoric_stress_tensor = xp->strength_data.deviatoric_stress_tensor_full;
}

/**
 * @brief Calculate time derivative of the deviatoric stress tensor.
 *
 * Hooke's law in rate form, dS/dt = 2 G (eps_dot - tr(eps_dot) I / 3), plus
 * the rotation terms that make it frame invariant.
 *
 * @param dS_dt The time derivative of the deviatoric stress tensor.
 * @param rotation_rate_tensor The rotation rate tensor.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void stress_tensor_compute_dS_dt(
    struct sym_matrix *dS_dt, const float rotation_rate_tensor[3][3],
    struct strength_state *s, const struct mat_params *mp) {

  const float shear_mod = mp->shear_mod;

  /* Convert strain rate to 3x3 float. */
  float strain_rate_tensor[3][3];
  get_matrix_from_sym_matrix(strain_rate_tensor, &s->strain_rate_tensor);

  /* Convert deviatoric stress to 3x3 float to compute the rotation term . */
  float deviatoric_stress_tensor[3][3];
  get_matrix_from_sym_matrix(deviatoric_stress_tensor,
                             &s->deviatoric_stress_tensor);

  /* Compute rotation term. */
  float rotation_term[3][3];
  strength_compute_rotation_term(rotation_term, rotation_rate_tensor,
                                 deviatoric_stress_tensor);

  /* Compute time derivative of the deviatoric stress tensor (Hooke's law). */
  float dS_dt_matrix[3][3];
  dS_dt_matrix[0][0] = 2.0f * shear_mod * strain_rate_tensor[0][0];
  dS_dt_matrix[0][1] = 2.0f * shear_mod * strain_rate_tensor[0][1];
  dS_dt_matrix[0][2] = 2.0f * shear_mod * strain_rate_tensor[0][2];
  dS_dt_matrix[1][0] = 2.0f * shear_mod * strain_rate_tensor[1][0];
  dS_dt_matrix[1][1] = 2.0f * shear_mod * strain_rate_tensor[1][1];
  dS_dt_matrix[1][2] = 2.0f * shear_mod * strain_rate_tensor[1][2];
  dS_dt_matrix[2][0] = 2.0f * shear_mod * strain_rate_tensor[2][0];
  dS_dt_matrix[2][1] = 2.0f * shear_mod * strain_rate_tensor[2][1];
  dS_dt_matrix[2][2] = 2.0f * shear_mod * strain_rate_tensor[2][2];

  /* Subtract trace from diagonal elements. */
  const float trace = strain_rate_tensor[0][0] + strain_rate_tensor[1][1] + strain_rate_tensor[2][2];
  dS_dt_matrix[0][0] -= 2.0f * shear_mod * trace / 3.f;
  dS_dt_matrix[1][1] -= 2.0f * shear_mod * trace / 3.f;
  dS_dt_matrix[2][2] -= 2.0f * shear_mod * trace / 3.f;

  /* Add rotation terms. */
  dS_dt_matrix[0][0] += rotation_term[0][0];
  dS_dt_matrix[0][1] += rotation_term[0][1];
  dS_dt_matrix[0][2] += rotation_term[0][2];
  dS_dt_matrix[1][0] += rotation_term[1][0];
  dS_dt_matrix[1][1] += rotation_term[1][1];
  dS_dt_matrix[1][2] += rotation_term[1][2];
  dS_dt_matrix[2][0] += rotation_term[2][0];
  dS_dt_matrix[2][1] += rotation_term[2][1];
  dS_dt_matrix[2][2] += rotation_term[2][2];

  /* Convert to sym_matrix. */
  get_sym_matrix_from_matrix(dS_dt, dS_dt_matrix);
}

/**
 * @brief Evolve the deviatoric stress tensor.
 *
 * @param deviatoric_stress_tensor The deviatoric stress tensor to be evolved.
 * @param dS_dt The time derivative of the deviatoric stress tensor.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 * @param dt_therm The time-step.
 */
__attribute__((always_inline)) INLINE static void stress_tensor_evolve_deviatoric_stress_tensor(
    struct sym_matrix *deviatoric_stress_tensor, const struct sym_matrix dS_dt,
    struct strength_state *s, const struct mat_params *mp,
    const float dt_therm) {

  /* Set to sym_matrix with all elements 0.f if the material is not solid. */
  if (strength_state_get_phase(s, mp) != mat_phase_solid) {
    zero_sym_matrix(deviatoric_stress_tensor);
    return;
  }

  /* Evolve deviatoric stress. */
  deviatoric_stress_tensor->xx += dS_dt.xx * dt_therm;
  deviatoric_stress_tensor->yy += dS_dt.yy * dt_therm;
  deviatoric_stress_tensor->zz += dS_dt.zz * dt_therm;
  deviatoric_stress_tensor->xy += dS_dt.xy * dt_therm;
  deviatoric_stress_tensor->xz += dS_dt.xz * dt_therm;
  deviatoric_stress_tensor->yz += dS_dt.yz * dt_therm;
}

/**
 * @brief Initialises the stress tensor  properties for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions or assignments between the particle
 * and extended particle fields.
 *
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static void strength_first_init_part_stress_tensor(
    struct part *restrict p, struct xpart *restrict xp) {

  zero_sym_matrix(&p->strength_data.dS_dt);
  zero_sym_matrix(&p->strength_data.deviatoric_stress_tensor);
  zero_sym_matrix(&xp->strength_data.deviatoric_stress_tensor_full);
  p->strength_data.is_yielding = 0;
}

#endif /* SWIFT_STRENGTH_STRESS_TENSOR_H */
