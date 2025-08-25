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
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Get tensile damage of particle at last drift time.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static float damage_get_tensile_damage(const struct part *restrict p) {

  return p->strength_data.tensile_damage;
}

/**
 * @brief Get tensile damage of particle at last kick time.
 *
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static float damage_get_tensile_damage_full(const struct xpart *restrict xp) {

  return xp->strength_data.tensile_damage_full;
}

/**
 * @brief Set drift-time tensile damage of particle.
 *
 * @param p The particle of interest.
 */
__attribute__((always_inline)) INLINE static void damage_set_tensile_damage(struct part *restrict p, const float tensile_damage) {

  p->strength_data.tensile_damage = tensile_damage;
}

/**
 * @brief Set kick-time tensile damage of particle.
 *
 * @param xp The extended data of the particle of interest.
 */
__attribute__((always_inline)) INLINE static void damage_set_tensile_damage_full(struct xpart *restrict xp, const float tensile_damage_full) {

  xp->strength_data.tensile_damage_full = tensile_damage_full;
}

/**
 * @brief Calculates the rate of cbrt(damage) accumulation due to tension.
 *
 * Benz&Asphaug95 tensile fracture model. Each particle is initialised with a
 * number of flaws, each with a corresponding activation threshold. A flaw is
 * activated if the local strain exceeds its activation threshold. Tensile
 * damage accumulates in a particle based on the number of curently-active flaws
 *  and an estimate of the crack velocity.
 *
 * Method parameters needed in material parameter file:
 * Strength:
 *     shear_mod: Shear modulus (Pa).
 *     bulk_mod: Bulk modulus (Pa).
 *
 * @param tensile_cbrtD_dt The rate of tensile cbrt(damage) accumulation.
 * @param number_of_activated_flaws The number of currently-active flaws.
 * @param number_of_flaws The total number of flaws.
 * @param activation_thresholds The activation thresholds fo flaws.
 * @param stress_tensor The stress tensor.
 * @param mat_id The material ID.
 * @param mass The particle mass.
 * @param density The density.
 * @param damage The damage.
 */
__attribute__((always_inline)) INLINE static void damage_tensile_compute_cbrtD_dt(
    float *tensile_cbrtD_dt, int *number_of_activated_flaws,  const int number_of_flaws, const float activation_thresholds[40], // ### Change this length
    const struct sym_matrix stress_tensor, const int mat_id, const float mass, const float density, const float damage) {

  *tensile_cbrtD_dt = 0.f;

  /* Tensile damage will only accumulate if a particle has flaws. */
  if (number_of_flaws == 0) {
    return;
  }

  /* Calculate maximum principal stress. */
  float principal_stress_eigen[3];
  sym_matrix_compute_eigenvalues(principal_stress_eigen, stress_tensor);
  float max_principal_stress = principal_stress_eigen[0];
  if (principal_stress_eigen[1] > max_principal_stress)
    max_principal_stress = principal_stress_eigen[1];
  if (principal_stress_eigen[2] > max_principal_stress)
    max_principal_stress = principal_stress_eigen[2];

  /* Tensile damage will only accumulate if particle is in tension. */
  if (max_principal_stress <= 0.f) {
    return;
  }

  /* Compute Young's modulus, E. */
  const float shear_mod = material_shear_mod(mat_id);
  const float bulk_mod = material_bulk_mod(mat_id);
  const float E = 9.f * bulk_mod * shear_mod / (3.f * bulk_mod + shear_mod);

  /* Compute local scalar strain, Benz&Asphaug95 Eqn. 21. */
  const float local_scalar_strain = max_principal_stress / ((1.f - damage) * E);

  /* The number of activated flaws is defined as the number of flaws for which
   * the local scalar strain has reached or exceeded the flaw's activation
   * threshold. */
  *number_of_activated_flaws = 0;
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
 * @brief Calculates the rate of damage accumulation due to tension
 *
 * @param tensile_dD_dt The rate of tensile damage accumulation.
 * @param p The particle of interest.
 * @param stress_tensor The stress tensor.
 * @param mat_id The material ID.
 * @param mass The particle mass.
 * @param density The density.
 * @param damage The damage.
 * @param tensile_damage The tensile damage.
 */
__attribute__((always_inline)) INLINE static void damage_tensile_compute_dD_dt(
    float *tensile_dD_dt, struct part *restrict p, const struct sym_matrix stress_tensor, const int mat_id, const float mass, const float density, const float damage, const float tensile_damage) {

  /* Particle flaws and their thresholds. */
  const float number_of_flaws = p->strength_data.number_of_flaws;
  float activation_thresholds[40]; //###
  memcpy(activation_thresholds, p->strength_data.activation_thresholds, sizeof(activation_thresholds));

  /* Compute the rate of cbrt(damage) accumulation due to tension. */
  float tensile_cbrtD_dt = 0.f;
  int number_of_activated_flaws = 0;
  damage_tensile_compute_cbrtD_dt(&tensile_cbrtD_dt, &number_of_activated_flaws, number_of_flaws,
                               activation_thresholds, stress_tensor, mat_id, mass, density, damage);

  /* Tensile damage is limited by the fraction of activated to total flaws. */
  if (tensile_damage < (float)number_of_activated_flaws / (float)number_of_flaws) {
    /* Chain rule d(D^(1/3))/dt = d(D^(1/3))/dD * dD/dt. */
    *tensile_dD_dt = 3.f * powf(tensile_damage, 2.f / 3.f) * tensile_cbrtD_dt;
  }
}

/**
 * @brief Steps tensile damage by applying time-step to a tensile_cbrtD_dt.
 *
 * @param tensile_damage The tensile damage.
 * @param tensile_cbrtD_dt The rate of tensile cbrt(damage) accumulation.
 * @param number_of_activated_flaws The number of currently-active flaws.
 * @param number_of_flaws The total number of flaws.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static void damage_tensile_apply_timestep_to_tensile_damage(
    float *tensile_damage, const float tensile_cbrtD_dt,
    const int number_of_activated_flaws, const int number_of_flaws, const float dt_therm) {

  // ### Can this be simplified based on damage_tensile_compute_dD_dt?

  /* Apply time-step. */
  float Delta_cbrtD = tensile_cbrtD_dt * dt_therm;

  /* Tensile damage is limited by the fraction of activated to total flaws. */
  const float max_cbrtD =
      cbrtf((float)number_of_activated_flaws / (float)number_of_flaws);
  const float max_Delta_cbrtD =
      fmaxf(0.f, max_cbrtD - cbrtf(*tensile_damage));
  if (Delta_cbrtD > max_Delta_cbrtD){
    Delta_cbrtD = max_Delta_cbrtD;
  }

  /* Update tensile damage. */
  const float evolved_D_cbrt = cbrtf(*tensile_damage) + Delta_cbrtD;
  *tensile_damage = powf(evolved_D_cbrt, 3.f);
}

/**
 * @brief Evolves tensile damage.
 *
 * Carries out all calculations required to step tensile damage in time.
 *
 * @param tensile_damage The tensile damage.
 * @param p The particle of interest.
 * @param stress_tensor The stress tensor.
 * @param mat_id The material ID.
 * @param mass The particle mass.
 * @param density The density.
 * @param damage The damage.
 * @param dt_therm The time-step duration.
 */
__attribute__((always_inline)) INLINE static void damage_tensile_evolve(
    float *tensile_damage, struct part *restrict p,  const struct sym_matrix stress_tensor, const int mat_id,
    const float mass, const float density, const float damage, const float dt_therm) {

  /* Particle flaws and their thresholds. */
  const float number_of_flaws = p->strength_data.number_of_flaws;
  float activation_thresholds[40];
  memcpy(activation_thresholds, p->strength_data.activation_thresholds, sizeof(activation_thresholds));

  float tensile_cbrtD_dt;
  int number_of_activated_flaws = 0;

  /* Compute the rate of cbrt(damage) accumulation due to tension. */
  damage_tensile_compute_cbrtD_dt(&tensile_cbrtD_dt, &number_of_activated_flaws, number_of_flaws,
                               activation_thresholds, stress_tensor, mat_id, mass, density, damage);

  /* Update tensile damage. */
  damage_tensile_apply_timestep_to_tensile_damage(tensile_damage, tensile_cbrtD_dt, number_of_activated_flaws, number_of_flaws, dt_therm);
}

#endif /* SWIFT_DAMAGE_TENSILE_BA95_H */
