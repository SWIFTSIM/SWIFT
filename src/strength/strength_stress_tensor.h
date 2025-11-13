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
 * @brief Hooke's Law model for elastc stress
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"
#include "strength_utilities.h"

/**
 * @brief Computes the stress tensor time-step of a given particle.
 *
 * Calculates a time-step based on the particle's rate of elastic stress
 * accumulation. If this time-step is smaller than dt_cfl, dt_cfl gets
 * overwritten to this damage time-step.
 *
 * @param dt_cfl The hydro (+ strength) time-step.
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static void strength_compute_timestep_stress_tensor(
    float *dt_cfl, const struct part *restrict p, const struct hydro_props *restrict hydro_properties) {

  const float elastic_timestep_factor = hydro_properties->CFL_condition; // ### Set as same as CFL factor for now. Treat this similarly to CFL
  const float norm_dS_dt = norm_sym_matrix(&p->strength_data.dS_dt);
  const float shear_mod = material_shear_mod(p->mat_id);

  if (norm_dS_dt * *dt_cfl > elastic_timestep_factor * shear_mod) {
    *dt_cfl = elastic_timestep_factor * shear_mod / norm_dS_dt;
  }
}

/**
 * @brief Updates the max wave speed of solids with the longitudinal wave speed.
 *
 * @param wave_speed The wave speed to be updated.
 * @param p The particle of interest.
 * @param soundspeed The sound speed.
 * @param density The sound density.
 */
__attribute__((always_inline)) INLINE static void
strength_compute_max_wave_speed_stress_tensor(float *wave_speed, const struct part *restrict p, const float soundspeed, const float density) {
  if (p->phase_state == mat_phase_state_solid) {
    /* Speed of longitudinal elastic wave. */
    const float shear_mod = material_shear_mod(p->mat_id);
    *wave_speed = sqrtf(soundspeed * soundspeed + (4.f / 3.f) * shear_mod / density);
  }
}

/**
 * @brief Compute the stress tensor.
 *
 * Constructs the full stress tensor by combining the deviatoric stress
 * and hydrostatic pressure. Applies any method-specific modifications,
 * such due to the accumulation of damage.
 *
 * @param p The particle of interest.
 * @param pressure The pressure.
 */
__attribute__((always_inline)) INLINE static void strength_compute_stress_tensor(
    struct part *restrict p, const float pressure) {

  /* Construct stress tensor. */
  p->strength_data.stress_tensor = p->strength_data.deviatoric_stress_tensor;
  p->strength_data.stress_tensor.xx -= pressure;
  p->strength_data.stress_tensor.yy -= pressure;
  p->strength_data.stress_tensor.zz -= pressure;

  /* Apply modifications due to damage. */
  const float damage = strength_get_damage(p);
  const struct sym_matrix damaged_deviatoric_stress_tensor = yield_compute_damaged_deviatoric_stress_tensor(p->strength_data.deviatoric_stress_tensor, damage);
  damage_compute_stress_tensor(&p->strength_data.stress_tensor, damaged_deviatoric_stress_tensor, pressure, damage);
}

/**
 * @brief Calculates the pairwise stress tensors for the force interaction.
 *
 * The stress tensors used for the force interaction between a specific pair of
 * particles. These differ from the particle's own stress tensor, since they
 * factor in the phases of the two particles as well as the contribution of
 * artificial stress for the pairwie interaction.
 *
 * @param pairwise_stress_tensor_i Stress tensor of particle i for its interactiion with j.
 * @param pairwise_stress_tensor_j Stress tensor of particle j for its interactiion with i.
 * @param pi First particle.
 * @param pj Second particle.
 * @param r The particle separation.
 */
__attribute__((always_inline)) INLINE static void
strength_set_pairwise_stress_tensors(float pairwise_stress_tensor_i[3][3],
                                     float pairwise_stress_tensor_j[3][3],
                                     const struct part *restrict pi,
                                     const struct part *restrict pj,
                                     const float r) {

  /* Only overwrite the fluid pairwise stress tensors if both particles are solid. */
  if ((pi->phase_state == mat_phase_state_solid) &&
      (pj->phase_state == mat_phase_state_solid)) {

    /* Get stress tensors. */
    get_matrix_from_sym_matrix(pairwise_stress_tensor_i, &pi->strength_data.stress_tensor);
    get_matrix_from_sym_matrix(pairwise_stress_tensor_j, &pj->strength_data.stress_tensor);

    /* Apply artificial stress. */
    artif_stress_apply_artif_stress_to_pairwise_stress_tensors(pairwise_stress_tensor_i,
                                                               pairwise_stress_tensor_j, pi, pj, r);
  }
}


/**
 * @brief Sets the values of additional particle stress tensor properties at a
 * kick time
 *
 * @param p The particle of interest.
 * @param xp The extended data of this particle.
 */__attribute__((always_inline)) INLINE static void strength_reset_predicted_values_stress_tensor(
    struct part *restrict p, const struct xpart *restrict xp) {

  p->strength_data.deviatoric_stress_tensor = xp->strength_data.deviatoric_stress_tensor_full;
}

/**
 * @brief Calculate time derivative of the deviatoric stress tensor.
 *
 * ### Description of how this ties in with Hooke's law
 *
 * @param p The particle of interest.
 * @param dv The velocity gradient dv/dr.
 */
__attribute__((always_inline)) INLINE static void stress_tensor_compute_dS_dt(struct part *restrict p, const float dv[3][3]) {

  const float shear_mod = material_shear_mod(p->mat_id);
  float strain_rate_tensor[3][3], rotation_rate_tensor[3][3], rotation_term[3][3];

  /* Compute strain rate and rotation rate. */
  strength_compute_strain_rate_tensor(strain_rate_tensor, dv);
  strength_compute_rotation_rate_tensor(rotation_rate_tensor, dv);

  /* Convert deviatoric stress to 3x3 float to compute the rotation term . */
  float deviatoric_stress_tensor[3][3];
  get_matrix_from_sym_matrix(deviatoric_stress_tensor,
                             &p->strength_data.deviatoric_stress_tensor);

  /* Compute rotation term. */
  strength_compute_rotation_term(rotation_term, rotation_rate_tensor,
                                 deviatoric_stress_tensor);

  /* Compute time derivative of the deviatoric stress tensor (Hooke's law). */
  float dS_dt[3][3];
  dS_dt[0][0] = 2.0f * shear_mod * strain_rate_tensor[0][0] + rotation_term[0][0];
  dS_dt[0][1] = 2.0f * shear_mod * strain_rate_tensor[0][1] + rotation_term[0][1];
  dS_dt[0][2] = 2.0f * shear_mod * strain_rate_tensor[0][2] + rotation_term[0][2];
  dS_dt[1][0] = 2.0f * shear_mod * strain_rate_tensor[1][0] + rotation_term[1][0];
  dS_dt[1][1] = 2.0f * shear_mod * strain_rate_tensor[1][1] + rotation_term[1][1];
  dS_dt[1][2] = 2.0f * shear_mod * strain_rate_tensor[1][2] + rotation_term[1][2];
  dS_dt[2][0] = 2.0f * shear_mod * strain_rate_tensor[2][0] + rotation_term[2][0];
  dS_dt[2][1] = 2.0f * shear_mod * strain_rate_tensor[2][1] + rotation_term[2][1];
  dS_dt[2][2] = 2.0f * shear_mod * strain_rate_tensor[2][2] + rotation_term[2][2];

  dS_dt[0][0] -= 2.0f * shear_mod * (strain_rate_tensor[0][0] + strain_rate_tensor[1][1] +
                                     strain_rate_tensor[2][2]) / 3.f;
  dS_dt[1][1] -= 2.0f * shear_mod * (strain_rate_tensor[0][0] + strain_rate_tensor[1][1] +
                                     strain_rate_tensor[2][2]) / 3.f;
  dS_dt[2][2] -= 2.0f * shear_mod * (strain_rate_tensor[0][0] + strain_rate_tensor[1][1] +
                                     strain_rate_tensor[2][2]) / 3.f;

  /* Update sym_matrix particle property. */
  get_sym_matrix_from_matrix(&p->strength_data.dS_dt, dS_dt);
}

/**
 * @brief Evolve the deviatoric stress tensor.
 *
 * @param deviatoric_stress_tensor the deviatoric stress.
 * @param p The particle of interest.
 * @param phase_state The phase state.
 * @param dt_therm The time-step.
 */
__attribute__((always_inline)) INLINE static void stress_tensor_evolve_deviatoric_stress_tensor(
    struct sym_matrix *deviatoric_stress_tensor, struct part *restrict p, const int phase_state,
    float dt_therm) {

  /* Return sym_matrix with all elements 0.f if the material is not solid. */
  if (phase_state != mat_phase_state_solid) {
    return zero_sym_matrix(deviatoric_stress_tensor);
  }

  /* Evolve deviatoric stress. */
  deviatoric_stress_tensor->xx += p->strength_data.dS_dt.xx * dt_therm;
  deviatoric_stress_tensor->yy += p->strength_data.dS_dt.yy * dt_therm;
  deviatoric_stress_tensor->zz += p->strength_data.dS_dt.zz * dt_therm;
  deviatoric_stress_tensor->xy += p->strength_data.dS_dt.xy * dt_therm;
  deviatoric_stress_tensor->xz += p->strength_data.dS_dt.xz * dt_therm;
  deviatoric_stress_tensor->yz += p->strength_data.dS_dt.yz * dt_therm;
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

  zero_sym_matrix(&p->strength_data.deviatoric_stress_tensor);
  zero_sym_matrix(&xp->strength_data.deviatoric_stress_tensor_full);
}

#endif /* SWIFT_STRENGTH_STRESS_TENSOR_H */
