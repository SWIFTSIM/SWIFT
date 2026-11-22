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
#ifndef SWIFT_DAMAGE_TENSILE_BENZ_ASPHAUG_H
#define SWIFT_DAMAGE_TENSILE_BENZ_ASPHAUG_H

/**
 * @file strength/damage/damage_tensile/damage_tensile_benz_asphaug.h
 * Benz&Asphaug95 tensile fracture model.
 */

#include "const.h"
#include "math.h"

/**
 * @brief Sets the tensile damage properties of the strength state from the particle.
 *
 * @param s The strength state of the particle.
 * @param p The particle of interest.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void damage_tensile_set_state(
    struct strength_state *s, const struct part *restrict p,
    const struct mat_params *mp) {

  s->number_of_flaws = p->strength_data.number_of_flaws;
  const struct strength_flaws *all_flaws = strength_get_flaws();
  /* View into all_flaws->activation_thresholds, so that this particle's flaws
   * are s->activation_thresholds[0, ..., s->number_of_flaws - 1]. */
  s->activation_thresholds =
      &all_flaws->activation_thresholds[p->strength_data.first_flaw_index];
}

/**
 * @brief Checks the flaws of a particle read from the ICs.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static void damage_tensile_first_init_part(
    const struct part *restrict p) {

  const long long first = p->strength_data.first_flaw_index;
  const int number = p->strength_data.number_of_flaws;
  const long long total = strength_get_flaws()->total_number_of_flaws;

  if ((number < 0) || (first < 0) || (first + number > total))
    error("Particle %lld has %d flaws starting at index %lld, but the ICs "
          "only have %lld flaws.",
          p->id, number, first, total);
}

/**
 * @brief Calculates the rate of cbrt(damage) accumulation due to tension.
 *
 * Benz&Asphaug95 tensile fracture model. Each particle is initialised with a
 * number of flaws, each with a corresponding activation threshold. A flaw is
 * activated if the local strain exceeds its activation threshold. Tensile
 * damage accumulates in a particle based on the number of currently-active flaws
 * and an estimate of the crack velocity.
 *
 * Method parameters needed in material parameter file:
 * Material:
 *     shear_mod: Shear modulus (Pa).
 *     bulk_mod: Bulk modulus (Pa).
 * DamageTensileBenzAsphaug:
 *     crack_velocity_factor: Crack velocity as a fraction of the longitudinal
 *     wave speed (default 0.4).
 *
 * @param damage_rates The rates of damage accumulation.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void damage_tensile_compute_cbrtD_dt(
    struct strength_damage_rates *damage_rates, struct strength_state *s,
    const struct mat_params *mp) {

  /* Particle properties. */
  const float mass = s->mass;
  const float density = s->density;

  /* Particle flaws and their thresholds. */
  const int number_of_flaws = s->number_of_flaws;
  const float *activation_thresholds = s->activation_thresholds;

  /* Tensile damage will only accumulate if a particle has flaws. */
  if (number_of_flaws == 0) {
    return;
  }

  /* If damage is already 1, no accumulation. */
  if (s->damage.total >= 1.f || s->damage.tensile >= 1.f) {
    return;
  }

  /* Calculate maximum principal stress. */
  const float *principal_stress_eigen = strength_state_get_principal_stresses(s, mp);
  const float max_principal_stress = fmaxf(principal_stress_eigen[0],
                                      fmaxf(principal_stress_eigen[1], principal_stress_eigen[2]));

  /* Tensile damage will only accumulate if particle is in tension. */
  if (max_principal_stress <= 0.f) {
    return;
  }

  /* Compute Young's modulus, E. */
  const float shear_mod = mp->shear_mod;
  const float bulk_mod = mp->bulk_mod;
  const float E = 9.f * bulk_mod * shear_mod / (3.f * bulk_mod + shear_mod);

  /* Compute local scalar strain, Benz&Asphaug95 Eqn. 21. Damage reduces the
   * stress for a given strain by the damage factor, so divide it back out. */
  const float local_scalar_strain =
      max_principal_stress / (yield_damage_factor(s, mp) * E);

  /* The number of activated flaws is defined as the number of flaws for which
   * the local scalar strain has reached or exceeded the flaw's activation
   * threshold. */
  int number_of_activated_flaws = 0;
  for (int i = 0; i < number_of_flaws; i++) {
    if (local_scalar_strain > activation_thresholds[i]) {
      number_of_activated_flaws += 1;
    }
  }

  /* Tensile damage is limited by the fraction of activated to total flaws. */
  damage_rates->tensile_max = (float)number_of_activated_flaws / (float)number_of_flaws;

  /* No accumulation once the damage has reached that limit. */
  if (s->damage.tensile >= damage_rates->tensile_max) {
    return;
  }

  /* Estimate the crack velocity based on the longitudinal wave speed,
   * Schäfer+2016 Eqn. 24. */
  const float longitudinal_wave_speed = sqrtf((bulk_mod + (4.f / 3.f) * (1.f - s->damage.total) * shear_mod) / density);
  const float crack_velocity = mp->damage_tensile_benz_asphaug_crack_velocity_factor * longitudinal_wave_speed;

  /* Compute the rate of cbrt(damage) accumulation due to tension,
   * Benz&Asphaug95 Eqn. 15. */
  damage_rates->tensile_cbrtD_dt = (float)number_of_activated_flaws * crack_velocity * cbrtf(density / mass);
}

/**
 * @brief Evolves tensile damage.
 *
 * @param d The damage of the particle.
 * @param damage_rates The rates of damage accumulation.
 * @param dt_therm The time-step duration.
 *
 * @return The evolved tensile damage.
 */
__attribute__((always_inline)) INLINE static float damage_tensile_evolve(
    struct strength_damage *d, const struct strength_damage_rates *damage_rates,
    const float dt_therm) {

  /* Tensile damage only accumulates with a positive rate. */
  if (damage_rates->tensile_cbrtD_dt <= 0.f) {
    return d->tensile;
  }

  /* Apply time-step. */
  float Delta_cbrtD = damage_rates->tensile_cbrtD_dt * dt_therm;

  /* Tensile damage is limited by the fraction of activated to total flaws. */
  /* See B&A99 and also Schafer16 for correction of cbrt placement here */
  const float max_Delta_cbrtD =
      cbrtf(damage_rates->tensile_max) - cbrtf(d->tensile);

  /* Already at or above the limit. */
  if (max_Delta_cbrtD <= 0.f) {
    return d->tensile;
  }

  /* Reaching the limit, so set damage to limit exactly. */
  if (Delta_cbrtD >= max_Delta_cbrtD) {
    d->tensile = damage_rates->tensile_max;
    return d->tensile;
  }

  /* Update tensile damage. */
  const float evolved_cbrtD = cbrtf(d->tensile) + Delta_cbrtD;
  d->tensile = fminf(powf(evolved_cbrtD, 3.f), 1.f);

  return d->tensile;
}

/**
 * @brief Limits the time-step by the growth of cbrt(tensile damage).
 *
 * @param dt_cfl The hydro + strength time-step.
 * @param damage_rates The rates of damage accumulation.
 * @param strength_props The material strength parameters.
 */
__attribute__((always_inline)) INLINE static void damage_tensile_compute_timestep(
    float *dt_cfl, const struct strength_damage_rates *damage_rates,
    const struct strength_props *strength_props) {

  /* Time for cbrt(tensile damage) to grow by damage_tensile_timestep_factor. */
  if (damage_rates->tensile_cbrtD_dt > 0.f) {
    *dt_cfl = fminf(*dt_cfl, strength_props->damage_tensile_timestep_factor /
                                 damage_rates->tensile_cbrtD_dt);
  }
}

#endif /* SWIFT_DAMAGE_TENSILE_BENZ_ASPHAUG_H */
