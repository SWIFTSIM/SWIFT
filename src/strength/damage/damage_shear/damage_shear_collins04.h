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
#ifndef SWIFT_DAMAGE_SHEAR_COLLINS04_H
#define SWIFT_DAMAGE_SHEAR_COLLINS04_H

/**
 * @file strength/damage/damage_shear/damage_shear_collins04.h
 * Collins+2004 shear fracture model.
 */

#include "const.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Get shear damage of particle at last drift time.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float damage_get_shear_damage(const struct part *restrict p) {

  return p->strength_data.shear_damage;
}

/**
 * @brief Get shear damage of particle at last kick time.
 *
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static float damage_get_shear_damage_full(const struct xpart *restrict xp) {

  return xp->strength_data.shear_damage_full;
}

/**
 * @brief Set drift-time shear damage of particle.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static void damage_set_shear_damage(struct part *restrict p, const float shear_damage) {

  p->strength_data.shear_damage = shear_damage;
}

/**
 * @brief Set kick-time shear damage of particle.
 *
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static void damage_set_shear_damage_full(struct xpart *restrict xp, const float shear_damage_full) {

  xp->strength_data.shear_damage_full = shear_damage_full;
}

/**
 * @brief Compute the plastic strain at failure.
 *
 * Method parameters needed in material parameter file:
 * DamageShearCollins:
 *     brittle_to_ductile_pressure: brittle--ductile transition pressure (Pa).
 *     brittle_to_plastic_pressure: brittle--plastic transition pressure (Pa).
 *
 * @param mat_id The material ID.
 * @param pressure The pressure.
 */
__attribute__((always_inline)) INLINE static float
damage_shear_compute_plastic_strain_at_failure(const int mat_id, const float pressure) {

  /* Method parameters. */
  const float brittle_to_ductile_pressure =
      material_brittle_to_ductile_pressure(mat_id);
  const float brittle_to_plastic_pressure =
      material_brittle_to_plastic_pressure(mat_id);

  float slope, intercept;
  if (pressure <= 0.f) {
    slope = 0.f;
    intercept = 0.01f;

  } else if (pressure < brittle_to_ductile_pressure) {
    /* Transitions between 0.01 and 0.05 in brittle regime. */
    slope = (0.05f - 0.01f) / brittle_to_ductile_pressure;
    intercept = 0.01f;

  } else if (pressure < brittle_to_plastic_pressure) {
    /* Transitions between 0.05 and 0.1 in semi-brittle regime. */
    slope = (0.1f - 0.05f) /
            (brittle_to_plastic_pressure - brittle_to_ductile_pressure);
    intercept = 0.05f - slope * brittle_to_ductile_pressure;

  } else {
    /* Transitions between 0.1 and 1 in plastic regime. */
    slope = (1.f - 0.1f) / brittle_to_plastic_pressure;
    intercept = 0.1f - slope * brittle_to_plastic_pressure;
  }

  return slope * pressure + intercept;
}

/**
 * @brief Compute the rate of shear damage accumulation.
 *
 * @param shear_dD_dt The rate of shear damage accumulation.
 * @param is_above_yield_criterion Whether the yield criterion is reached.
 * @param strain_rate_tensor The strain rate tensor.
 * @param mat_id The material ID.
 * @param pressure The pressure.
 * @param shear_damage The shear damage.
 */
__attribute__((always_inline)) INLINE static void damage_shear_compute_dD_dt(
    float *shear_dD_dt, const int is_above_yield_criterion, const struct sym_matrix strain_rate_tensor,
    const int mat_id, const float pressure, const float shear_damage) {

  *shear_dD_dt = 0.f;

  /* If shear damage is already 1, no accumulation. */
  if (shear_damage >= 1.f) {
    return;
  }

  /* Shear damage not accumulated if yield criterion is not reached. */
  if (!is_above_yield_criterion) {
    return;
  }

  /* Calculate invariant of strain rate tensor. While yielding, assume that the
   * plastic strain rate is approx. the total strain rate. */
  const float E_dot = sqrtf(strength_compute_sym_matrix_J_2(strain_rate_tensor));

  /* Calculate the plastic strain at failure. */
  const float plastic_strain_at_failure =
      damage_shear_compute_plastic_strain_at_failure(mat_id, pressure);

  /* Calculate the rate of shear damage accumulation. */
  *shear_dD_dt = E_dot / plastic_strain_at_failure;
}

/**
 * @brief Steps shear damage by applying time-step to a shear_dD_dt.
 *
 * @param shear_damage The shear damage.
 * @param shear_dD_dt The rate of shear damage accumulation.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static void damage_shear_apply_timestep_to_shear_damage(
    float *shear_damage, const float shear_dD_dt, const float dt_therm) {

  if (shear_dD_dt <= 0.f) {
    return;
  }

  /* Update shear damage. */
  *shear_damage = fminf(*shear_damage + shear_dD_dt * dt_therm, 1.f);
}

/**
 * @brief Evolves shear damage.
 *
 * Carries out all calculations required to step shear damage in time.
 *
 * @param shear_damage The shear damage.
 * @param is_above_yield_criterion Whether the yield criterion is reached.
 * @param strain_rate_tensor The strain rate tensor.
 * @param mat_id The material ID.
 * @param pressure The pressure.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static void damage_shear_evolve(
    float *shear_damage, const int is_above_yield_criterion, const struct sym_matrix strain_rate_tensor,
    const int mat_id, const float pressure, const float dt_therm) {

  float shear_dD_dt;

  /* Compute the rate of shear damage accumulation. */
  damage_shear_compute_dD_dt(&shear_dD_dt, is_above_yield_criterion, strain_rate_tensor, mat_id, pressure, *shear_damage);

  /* Update shear damage. */
  damage_shear_apply_timestep_to_shear_damage(shear_damage, shear_dD_dt, dt_therm);
}

#endif /* SWIFT_DAMAGE_SHEAR_COLLINS04_H */
