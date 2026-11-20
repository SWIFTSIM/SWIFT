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
    const struct mat_params *mp) {

  /* Whether the yield stress limited S in the last drift. Only set in the
   * drift, like the damage rates, which are computed from the drifted values. */
  s->is_yielding = p->strength_data.is_yielding;
}

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
 * @param damage_rates The rates of damage accumulation.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void damage_shear_compute_dD_dt(
    struct strength_damage_rates *damage_rates, struct strength_state *s,
    const struct mat_params *mp) {

  /* If damage is already 1, no accumulation. */
  if (s->damage.total >= 1.f || s->damage.shear >= 1.f) {
    return;
  }

  /* Shear damage only accumulates if the yield stress limited the deviatoric
   * stress tensor in the last drift. */
  if (!s->is_yielding) {
    return;
  }

  /* Strain rate invariant, Collins+2004 Eqn. A7. While yielding, assume that
   * all of the deviatoric strain rate is plastic. */
  const float E_dot = sqrtf(strength_state_get_strain_rate_J_2(s, mp));

  /* Calculate the plastic strain at failure. */
  const float plastic_strain_at_failure =
      damage_shear_compute_plastic_strain_at_failure(s, mp);

  /* Calculate the rate of shear damage accumulation. */
  damage_rates->shear_dD_dt = E_dot / plastic_strain_at_failure;
}

/**
 * @brief Evolves shear damage.
 *
 * @param d The damage of the particle.
 * @param damage_rates The rates of damage accumulation.
 * @param dt_therm The time-step duration.
 *
 * @return The evolved shear damage.
 */
__attribute__((always_inline)) INLINE static float damage_shear_evolve(
    struct strength_damage *d, const struct strength_damage_rates *damage_rates,
    const float dt_therm) {

  if (damage_rates->shear_dD_dt <= 0.f) {
    return d->shear;
  }

  /* Update shear damage. */
  d->shear = fminf(fmaxf(d->shear + damage_rates->shear_dD_dt * dt_therm, 0.f), 1.f);

  return d->shear;
}

/**
 * @brief Limits the time-step by the growth of shear damage.
 *
 * @param dt_cfl The hydro + strength time-step.
 * @param damage_rates The rates of damage accumulation.
 * @param strength_props The material strength parameters.
 */
__attribute__((always_inline)) INLINE static void damage_shear_compute_timestep(
    float *dt_cfl, const struct strength_damage_rates *damage_rates,
    const struct strength_props *strength_props) {

  /* Time for shear damage to grow by damage_shear_timestep_factor. */
  if (damage_rates->shear_dD_dt > 0.f) {
    *dt_cfl = fminf(*dt_cfl, strength_props->damage_shear_timestep_factor /
                                 damage_rates->shear_dD_dt);
  }
}

#endif /* SWIFT_DAMAGE_SHEAR_COLLINS_H */
