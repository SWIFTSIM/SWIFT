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
#ifndef SWIFT_DAMAGE_SHEAR_COLLINS_H
#define SWIFT_DAMAGE_SHEAR_COLLINS_H

/**
 * @file strength/damage/damage_shear/damage_shear_collins.h
 * Collins+2004 shear fracture model.
 */

#include "const.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Sets the shear damage properties of the strength state from the particle.
 *
 * @param s The strength state of the particle.
 * @param p The particle of interest.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void damage_shear_set_state(
    struct strength_state *s, const struct part *restrict p,
    const struct mat_params *mp) {}

/**
 * @brief Compute the plastic strain at failure.
 *
 * Method parameters needed in material parameter file:
 * DamageShearCollins:
 *     brittle_to_ductile_pressure: brittle--ductile transition pressure (Pa).
 *     brittle_to_plastic_pressure: brittle--plastic transition pressure (Pa).
 *     plastic_strain_at_failure: Plastic strain at failure at zero pressure,
 *     at the two transition pressures, and at twice the brittle--plastic
 *     transition pressure (default [0.01, 0.05, 0.1, 1]).
 *
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static float
damage_shear_compute_plastic_strain_at_failure(struct strength_state *s,
                                               const struct mat_params *mp) {

  /* Method parameters. */
  const float brittle_to_ductile_pressure =
      mp->damage_shear_collins_brittle_to_ductile_pressure;
  const float brittle_to_plastic_pressure =
      mp->damage_shear_collins_brittle_to_plastic_pressure;
  const float *plastic_strain_at_failure = mp->damage_shear_collins_plastic_strain_at_failure;
  const float pressure = strength_state_get_pressure(s, mp);

  float slope, intercept;
  if (pressure <= 0.f) {
    slope = 0.f;
    intercept = plastic_strain_at_failure[0];

  } else if (pressure < brittle_to_ductile_pressure) {
    /* Transitions between plastic_strain_at_failure[0] and plastic_strain_at_failure[1] in brittle regime. */
    slope = (plastic_strain_at_failure[1] - plastic_strain_at_failure[0]) / brittle_to_ductile_pressure;
    intercept = plastic_strain_at_failure[0];

  } else if (pressure < brittle_to_plastic_pressure) {
    /* Transitions between plastic_strain_at_failure[1] and plastic_strain_at_failure[2] in semi-brittle regime. */
    slope = (plastic_strain_at_failure[2] - plastic_strain_at_failure[1]) /
            (brittle_to_plastic_pressure - brittle_to_ductile_pressure);
    intercept = plastic_strain_at_failure[1] - slope * brittle_to_ductile_pressure;

  } else {
    /* Transitions between plastic_strain_at_failure[2] and plastic_strain_at_failure[3] in plastic regime. */
    slope = (plastic_strain_at_failure[3] - plastic_strain_at_failure[2]) / brittle_to_plastic_pressure;
    intercept = plastic_strain_at_failure[2] - slope * brittle_to_plastic_pressure;
  }

  return slope * pressure + intercept;
}

/**
 * @brief Compute the rate of shear damage accumulation.
 *
 * @param shear_dD_dt The rate of shear damage accumulation.
 * @param shear_damage The shear damage.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void damage_shear_compute_dD_dt(
    float *shear_dD_dt, const float shear_damage,
    struct strength_state *s, const struct mat_params *mp) {

  *shear_dD_dt = 0.f;

  /* If shear damage is already 1, no accumulation. */
  if (shear_damage >= 1.f) {
    return;
  }

  /* Shear damage not accumulated if yield criterion is not reached. */
  if (!yield_check_yield_criterion(s, mp)) {
    return;
  }

  /* Calculate invariant of strain rate tensor. While yielding, assume that the
   * plastic strain rate is approx. the total strain rate. */
  const float E_dot = sqrtf(strength_state_get_strain_rate_J_2(s, mp));

  /* Calculate the plastic strain at failure. */
  const float plastic_strain_at_failure =
      damage_shear_compute_plastic_strain_at_failure(s, mp);

  /* Calculate the rate of shear damage accumulation. */
  *shear_dD_dt = E_dot / plastic_strain_at_failure;
}

/**
 * @brief Steps shear damage by applying time-step to a shear_dD_dt.
 *
 * @param shear_damage The shear damage.
 * @param shear_dD_dt The rate of shear damage accumulation.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static void damage_shear_apply_timestep_to_shear_damage(
    float *shear_damage, const float shear_dD_dt,
    struct strength_state *s, const struct mat_params *mp,
    const float dt_therm) {

  if (shear_dD_dt <= 0.f) {
    return;
  }

  /* Update shear damage. */
  *shear_damage = fminf(fmaxf(*shear_damage + shear_dD_dt * dt_therm, 0.f), 1.f);
}

/**
 * @brief Evolves shear damage.
 *
 * Carries out all calculations required to step shear damage in time.
 *
 * @param d The damage of the particle. Only the shear damage is changed.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 * @param dt_therm The time-step duration.
 *
 * @return The evolved shear damage.
 */
__attribute__((always_inline)) INLINE static float damage_shear_evolve(
    struct strength_damage *d, struct strength_state *s,
    const struct mat_params *mp, const float dt_therm) {

  float shear_dD_dt;

  /* Compute the rate of shear damage accumulation. */
  damage_shear_compute_dD_dt(&shear_dD_dt, d->shear, s, mp);

  /* Update shear damage. */
  damage_shear_apply_timestep_to_shear_damage(&d->shear, shear_dD_dt, s, mp, dt_therm);

  return d->shear;
}

/**
 * @brief Computes the rate of shear damage accumulation used to set the
 * damage time-step.
 *
 * Shear damage evolves D directly, so this is just its rate.
 *
 * @param d The damage of the particle.
 * @param damage_timestep_factor The growth in damage allowed per time-step.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 *
 * @return The rate of shear damage accumulation.
 */
__attribute__((always_inline)) INLINE static float damage_shear_compute_timestep_dD_dt(
    const struct strength_damage *d, const float damage_timestep_factor,
    struct strength_state *s, const struct mat_params *mp) {

  float shear_dD_dt;
  damage_shear_compute_dD_dt(&shear_dD_dt, d->shear, s, mp);

  return shear_dD_dt;
}

#endif /* SWIFT_DAMAGE_SHEAR_COLLINS_H */
