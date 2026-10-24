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
#ifndef SWIFT_DAMAGE_TENSILE_BA95_H
#define SWIFT_DAMAGE_TENSILE_BA95_H

/**
 * @file strength/damage/damage_tensile/damage_tensile_ba95.h
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
  s->activation_thresholds = p->strength_data.activation_thresholds;
}

/**
 * @brief Calculates the rate of cbrt(damage) accumulation due to tension.
 *
 * Benz&Asphaug95 tensile fracture model. Each particle is initialised with a
 * number of flaws, each with a corresponding activation threshold. A flaw is
 * activated if the local strain exceeds its activation threshold. Tensile
 * damage accumulates in a particle based on the number of curently-active flaws
 * and an estimate of the crack velocity.
 *
 * Method parameters needed in material parameter file:
 * Strength:
 *     shear_mod: Shear modulus (Pa).
 *     bulk_mod: Bulk modulus (Pa).
 *
 * @param tensile_cbrtD_dt The rate of tensile cbrt(damage) accumulation.
 * @param number_of_activated_flaws The number of currently-active flaws.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void damage_tensile_compute_cbrtD_dt(
    float *tensile_cbrtD_dt, int *number_of_activated_flaws,
    const struct strength_state *s, const struct mat_params *mp) {

  *tensile_cbrtD_dt = 0.f;
  *number_of_activated_flaws = 0;

  /* Particle properties. */
  const struct sym_matrix stress_tensor = s->stress_tensor;
  const float mass = s->mass;
  const float density = s->density;
  const float damage = s->damage.total;

  /* Particle flaws and their thresholds. */
  const int number_of_flaws = s->number_of_flaws;
  const float *activation_thresholds = s->activation_thresholds;

  /* Tensile damage will only accumulate if a particle has flaws. */
  if (number_of_flaws == 0) {
    return;
  }

  /* If damage is already 1, no accumulation. */
  if (damage == 1.f) {
    return;
  }

  /* Calculate maximum principal stress. */
  float principal_stress_eigen[3];
  sym_matrix_compute_eigenvalues(principal_stress_eigen, stress_tensor);
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

  /* Compute local scalar strain, Benz&Asphaug95 Eqn. 21. */
  const float local_scalar_strain = max_principal_stress / ((1.f - damage) * E);

  /* The number of activated flaws is defined as the number of flaws for which
   * the local scalar strain has reached or exceeded the flaw's activation
   * threshold. */
  // ### this can be done nicer if we sort thresholds in order and then break in this loop
  for (int i = 0; i < number_of_flaws; i++) {
    if (local_scalar_strain > activation_thresholds[i]) {
      *number_of_activated_flaws += 1;
    }
  }

  /* Estimate the crack velocity based on the longitudinal wave speed,
   * Schäfer+2016 Eqn. 24. */
  const float longitudinal_wave_speed = sqrtf((bulk_mod + (4.f / 3.f) * (1.f - damage) * shear_mod) / density);
  const float crack_velocity = 0.4f * longitudinal_wave_speed;

  /* Compute the rate of cbrt(damage) accumulation due to tension,
   * Benz&Asphaug95 Eqn. 15. */
  *tensile_cbrtD_dt = (float)*number_of_activated_flaws * crack_velocity * cbrtf(density / mass);
}

/**
 * @brief Steps tensile damage by applying time-step to a tensile_cbrtD_dt.
 *
 * @param tensile_damage The tensile damage.
 * @param tensile_cbrtD_dt The rate of tensile cbrt(damage) accumulation.
 * @param number_of_activated_flaws The number of currently-active flaws.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static void damage_tensile_apply_timestep_to_tensile_damage(
    float *tensile_damage, const float tensile_cbrtD_dt,
    const int number_of_activated_flaws, const struct strength_state *s,
    const struct mat_params *mp, const float dt_therm) {

  /* Tensile damage only accumulates if the particle has flaws and a positive
   * rate. */
  const int number_of_flaws = s->number_of_flaws;
  if (number_of_flaws == 0 || tensile_cbrtD_dt <= 0.f) {
    return;
  }

  /* Apply time-step. */
  float Delta_cbrtD = tensile_cbrtD_dt * dt_therm;

  /* Tensile damage is limited by the fraction of activated to total flaws. */
  /* See B&A99 and also Schafer16 for correction of cbrt placement here */
  const float max_cbrtD =
      cbrtf((float)number_of_activated_flaws / (float)number_of_flaws);
  const float max_Delta_cbrtD =
      fmaxf(0.f, max_cbrtD - cbrtf(*tensile_damage));
  if (Delta_cbrtD > max_Delta_cbrtD){
    Delta_cbrtD = max_Delta_cbrtD;
  }

  /* Update tensile damage. */
  const float evolved_D_cbrt = fmaxf(cbrtf(*tensile_damage) + Delta_cbrtD, 0.f);
  *tensile_damage = fminf(powf(evolved_D_cbrt, 3.f), 1.f);
}

/**
 * @brief Evolves tensile damage.
 *
 * Carries out all calculations required to step tensile damage in time.
 *
 * @param d The damage of the particle. Only the tensile damage is changed.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 * @param dt_therm The time-step duration.
 *
 * @return The evolved tensile damage.
 */
__attribute__((always_inline)) INLINE static float damage_tensile_evolve(
    struct strength_damage *d, const struct strength_state *s,
    const struct mat_params *mp, const float dt_therm) {

  float tensile_cbrtD_dt;
  int number_of_activated_flaws = 0;

  /* Compute the rate of cbrt(damage) accumulation due to tension. */
  damage_tensile_compute_cbrtD_dt(&tensile_cbrtD_dt, &number_of_activated_flaws, s, mp);

  /* Update tensile damage. */
  damage_tensile_apply_timestep_to_tensile_damage(&d->tensile, tensile_cbrtD_dt, number_of_activated_flaws, s, mp, dt_therm);

  return d->tensile;
}

/**
 * @brief Computes the rate of tensile damage accumulation used to set the
 * damage time-step.
 *
 * @param d The damage of the particle.
 * @param damage_timestep_factor The growth in damage allowed per time-step.
 * @param s The strength state of the particle.
 * @param mp The material parameters.
 *
 * @return The rate of tensile damage accumulation.
 */
__attribute__((always_inline)) INLINE static float damage_tensile_compute_timestep_dD_dt(
    const struct strength_damage *d, const float damage_timestep_factor,
    const struct strength_state *s, const struct mat_params *mp) {

  const float tensile_damage = d->tensile;

  /* If tensile damage is already 1, no accumulation. */
  if (tensile_damage >= 1.f) {
    return 0.f;
  }

  float tensile_cbrtD_dt = 0.f;
  int number_of_activated_flaws = 0;
  damage_tensile_compute_cbrtD_dt(&tensile_cbrtD_dt, &number_of_activated_flaws,
                                  s, mp); // ### Should this use damage or tensile_damage?

  if (tensile_cbrtD_dt <= 0.f) {
    return 0.f;
  }

  const float Delta_cbrtD =
      cbrtf(tensile_damage + damage_timestep_factor) - cbrtf(tensile_damage);
  return damage_timestep_factor * tensile_cbrtD_dt / Delta_cbrtD;
}

#endif /* SWIFT_DAMAGE_TENSILE_BA95_H */
