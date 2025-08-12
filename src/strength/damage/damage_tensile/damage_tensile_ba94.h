/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Thomas Sandnes (thomas.d.sandnes@durham.ac.uk)
 *               2024 Jacob Kegerreis (jacob.kegerreis@durham.ac.uk)
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
#ifndef SWIFT_DAMAGE_TENSILE_BA94_H
#define SWIFT_DAMAGE_TENSILE_BA94_H

/**
 * @file strength/damage/damage_tensile/damage_tensile_ba94.h
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"

/**
 * @brief Calculates the rate of cbrt(damage) accumulated due to tension
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void calculate_tensile_cbrtD_dt(
    float *tensile_cbrtD_dt, int *number_of_activated_flaws,  const int number_of_flaws, const float activation_thresholds[40], // ### Change this length
    const struct sym_matrix stress_tensor, const int mat_id, const float mass, const float density, const float damage) {

  *tensile_cbrtD_dt = 0.f;  
  // Damage can only accumulate if in tension (this will have to change if we
  // add fracture under compression)

  // If there are no flaws, we can not accumulate damage
  if (number_of_flaws > 0) {
      
    float principal_stress_eigen[3];
    sym_matrix_compute_eigenvalues(principal_stress_eigen, stress_tensor); 

    // Find max eignenvalue of stress_tensor_sigma:
    float max_principal_stress = principal_stress_eigen[0];
    if (principal_stress_eigen[1] > max_principal_stress)
      max_principal_stress = principal_stress_eigen[1];
    if (principal_stress_eigen[2] > max_principal_stress)
      max_principal_stress = principal_stress_eigen[2];
     
    if (max_principal_stress > 0.f) {
      const float shear_mod = material_shear_mod(mat_id);
      const float bulk_mod = material_bulk_mod(mat_id);

      // tensile damage
      const float E = 9.f * bulk_mod * shear_mod / (3.f * bulk_mod + shear_mod);

      const float local_scalar_strain =
          max_principal_stress / ((1.f - damage) * E);

      *number_of_activated_flaws = 0;
      for (int i = 0; i < number_of_flaws; i++) {
        if (local_scalar_strain > activation_thresholds[i]) {
          *number_of_activated_flaws += 1;
        }
      }

      // Schafer
      const float longitudinal_wave_speed = sqrtf(
          (bulk_mod + (4.f / 3.f) * (1.f - damage) * shear_mod) / density);
      const float crack_velocity = 0.4f * longitudinal_wave_speed;

      *tensile_cbrtD_dt =
          (float)*number_of_activated_flaws * crack_velocity * cbrtf(density / mass);
    }
  }
}

/**
 * @brief Calculates the rate of damage accumulated due to tension
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void calculate_tensile_dD_dt(
    struct part *restrict p, float *tensile_dD_dt, const struct sym_matrix stress_tensor, const int mat_id, const float mass, const float density, const float damage, const float tensile_damage) {

  const float number_of_flaws = p->strength_data.number_of_flaws;
  float activation_thresholds[40];
  memcpy(activation_thresholds, p->strength_data.activation_thresholds, sizeof(activation_thresholds));

  float tensile_cbrtD_dt = 0.f;
  int number_of_activated_flaws = 0;
  calculate_tensile_cbrtD_dt(&tensile_cbrtD_dt, &number_of_activated_flaws, number_of_flaws, 
                               activation_thresholds, stress_tensor, mat_id, mass, density, damage);
  
  if (tensile_damage < (float)number_of_activated_flaws / (float)number_of_flaws) {
    // Chain rule d(D^(1/3))/dt  = d(D^(1/3))/dD * dD/dt
    *tensile_dD_dt = 3.f * powf(tensile_damage, 2.f / 3.f) * tensile_cbrtD_dt;
  }
}

/**
 * @brief Steps particle tensile damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void step_damage_tensile(
    float *tensile_damage, const float tensile_cbrtD_dt, 
    const int number_of_activated_flaws, const int number_of_flaws, const float dt_therm) {

  float Delta_cbrtD = tensile_cbrtD_dt * dt_therm;

  const float max_cbrtD =
      cbrtf((float)number_of_activated_flaws / (float)number_of_flaws);
  const float max_Delta_cbrtD =
      max(0.f, max_cbrtD - cbrtf(*tensile_damage));

  if (Delta_cbrtD > max_Delta_cbrtD) Delta_cbrtD = max_Delta_cbrtD;

  const float evolved_D_cbrt = cbrtf(*tensile_damage) + Delta_cbrtD;
  *tensile_damage = powf(evolved_D_cbrt, 3.f);
}

/**
 * @brief Evolves particle tensile damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void evolve_damage_tensile(
    struct part *restrict p, float *tensile_damage, const struct sym_matrix stress_tensor, const int mat_id, 
    const float mass, const float density, const float damage, const float dt_therm) {

  const float number_of_flaws = p->strength_data.number_of_flaws;
  float activation_thresholds[40];
  memcpy(activation_thresholds, p->strength_data.activation_thresholds, sizeof(activation_thresholds));
        
  float tensile_cbrtD_dt;
  int number_of_activated_flaws;
  calculate_tensile_cbrtD_dt(&tensile_cbrtD_dt, &number_of_activated_flaws, number_of_flaws, 
                               activation_thresholds, stress_tensor, mat_id, mass, density, damage);
  step_damage_tensile(tensile_damage, tensile_cbrtD_dt, number_of_activated_flaws, number_of_flaws, dt_therm);
}

#endif /* SWIFT_DAMAGE_TENSILE_BA94_H */