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
#ifndef SWIFT_DAMAGE_H
#define SWIFT_DAMAGE_H

/**
 * @file strength/damage/damage.h
 * @brief Core damage scheme that combines multiple fracture models.
 */

#include "const.h"
#include <float.h>
#include "math.h"

/**
 * @brief Get damage of particle at last drift time.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static struct strength_damage strength_get_damage(const struct part *restrict p) {

  return p->strength_data.damage;
}

/**
 * @brief Get damage of particle at last kick time.
 *
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static struct strength_damage strength_get_damage_full(const struct xpart *restrict xp) {

  return xp->strength_data.damage_full;
}

/**
 * @brief Sets the damage properties of the strength state from the particle.
 *
 * @param s The strength state of the particle.
 * @param p The particle of interest.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void damage_set_state(
    struct strength_state *s, const struct part *restrict p,
    const struct mat_params *mp) {

  damage_tensile_set_state(s, p, mp);
  damage_shear_set_state(s, p, mp);
}

/**
 * @brief Limits the time-step by the damage time-step.
 *
 * @param dt_cfl The hydro (+ strength) time-step.
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static void strength_compute_timestep_damage(
    float *dt_cfl, const struct part *restrict p) {

  *dt_cfl = fminf(*dt_cfl, p->strength_data.dt_damage);
}

/**
 * @brief Computes the damage-modified stress tensor.
 *
 * @param stress_tensor The stress tensor.
 * @param damaged_deviatoric_stress_tensor The damaged_deviatoric stress tensor, already modified by damage.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void damage_compute_stress_tensor(
    struct sym_matrix *stress_tensor, const struct sym_matrix damaged_deviatoric_stress_tensor,
    struct strength_state *s, const struct mat_params *mp) {

  /* How damage affects the deviatoric stress tensor depends on the yield stress
   * method used. */
  *stress_tensor = damaged_deviatoric_stress_tensor;

  /* Damage weakens negative pressures so that fully damaged material cannot be
   * in tension. */
  const float pressure = strength_state_get_pressure(s, mp);
  const float effective_pressure = pressure > 0.f ? pressure : (1.f - s->damage.total) * pressure;
  stress_tensor->xx -= effective_pressure;
  stress_tensor->yy -= effective_pressure;
  stress_tensor->zz -= effective_pressure;
}

/**
 * @brief Sets the values of additional particle damage properties at a
 * kick time
 *
 * @param p The particle of interest.
 * @param xp The extended data of this particle.
 */
__attribute__((always_inline)) INLINE static void strength_reset_predicted_values_damage(
    struct part *restrict p, const struct xpart *restrict xp) {

  p->strength_data.damage = xp->strength_data.damage_full;
}

/**
 * @brief Evolves particle damage.
 *
 * @param d The damage, to be evolved.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static void damage_evolve(
    struct strength_damage *d, struct strength_state *s,
    const struct mat_params *mp, const float dt_therm) {

  /* Only solid particles accumulate damage. */
  if (strength_state_get_phase(s, mp) != mat_phase_solid) {
    return;
  }

  /* Evolve each source of damage. */
  const float tensile_damage = damage_tensile_evolve(d, s, mp, dt_therm);
  const float shear_damage = damage_shear_evolve(d, s, mp, dt_therm);

  /* Combine sources of damage. */
  d->total = fminf(tensile_damage + shear_damage, 1.f);
}

/**
 * @brief Evolves particle damage in the drift
 *
 * @param p The particle of interest.
 * @param s The strength state of the particle, from its drifted values.
 * @param mp The material parameters.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static void damage_predict_evolve(
    struct part *restrict p, struct strength_state *s,
    const struct mat_params *mp, const float dt_therm) {

  damage_evolve(&p->strength_data.damage, s, mp, dt_therm);
}

/**
 * @brief Evolves particle damage in the kick
 *
 * @param xp The extended data of the particle of interest.
 * @param s The strength state of the particle, from its kicked values.
 * @param mp The material parameters.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static void damage_kick_evolve(
    struct xpart *restrict xp, struct strength_state *s,
    const struct mat_params *mp, const float dt_therm) {

  damage_evolve(&xp->strength_data.damage_full, s, mp, dt_therm);
}

/**
 * @brief Calculate the damage time-step.
 *
 * The time for damage to grow by damage_timestep_factor.
 *
 * @param p The particle of interest.
 * @param s The strength state of the particle, from its drifted values.
 * @param mp The material parameters.
 * @param strength_props The material strength parameters.
 */
__attribute__((always_inline)) INLINE static void damage_compute_timescale(
    struct part *restrict p, struct strength_state *s,
    const struct mat_params *mp,
    const struct strength_props *strength_props) {

  /* Growth in damage used to set the time-step. */
  const float damage_timestep_factor = strength_props->damage_timestep_factor;

  /* Damage at drift time. */
  const struct strength_damage *d = &p->strength_data.damage;

  p->strength_data.dt_damage = FLT_MAX;

  /* Only solid particles that are not fully damaged can accumulate damage. */
  if ((strength_state_get_phase(s, mp) != mat_phase_solid) || (d->total >= 1.f)) {
    return;
  }

  /* Combined rate of total damage from both sources. */
  const float dD_dt =
      damage_tensile_compute_timestep_dD_dt(d, damage_timestep_factor, s, mp) +
      damage_shear_compute_timestep_dD_dt(d, damage_timestep_factor, s, mp);

  /* Time for total damage to grow by damage_timestep_factor. */
  if (dD_dt > 0.f) {
    p->strength_data.dt_damage = damage_timestep_factor / dD_dt;
  }
}

/**
 * @brief Initialises the damage properties for the first time
 *
 * This function is called only once just after the ICs have been
 * read in to do some conversions or assignments between the particle
 * and extended particle fields.
 *
 * @param p The particle of interest.
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static void strength_first_init_part_damage(
    struct part *restrict p, struct xpart *restrict xp) {

  const struct strength_damage no_damage = {0};
  p->strength_data.damage = no_damage;
  xp->strength_data.damage_full = no_damage;

  damage_tensile_first_init_part(p);
}

#endif /* SWIFT_DAMAGE_H */
