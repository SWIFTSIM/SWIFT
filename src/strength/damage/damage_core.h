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
#ifndef SWIFT_DAMAGE_H
#define SWIFT_DAMAGE_H

/**
 * @file strength/damage/damage.h
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"
#include "strength.h"

__attribute__((always_inline)) INLINE static float strength_get_damage(const struct part *restrict p) {

  return p->strength_data.damage;
}

__attribute__((always_inline)) INLINE static float strength_get_damage_full(const struct xpart *restrict xp) {

  return xp->strength_data.damage_full;
}

__attribute__((always_inline)) INLINE static void strength_set_damage(struct part *restrict p, const float damage) {

  p->strength_data.damage = damage;
}

__attribute__((always_inline)) INLINE static void strength_set_damage_full(struct xpart *restrict xp, const float damage_full) {

  xp->strength_data.damage_full = damage_full;
}

__attribute__((always_inline)) INLINE static void strength_compute_timestep_damage(
    const struct part *restrict p, float *dt_cfl) {

  const float dD_dt = p->strength_data.dD_dt;

  const float damage_timestep_factor = 0.01f; // ### Hardcoded for now. Treat this similarly to CFL
    
  if (dD_dt * *dt_cfl > damage_timestep_factor) {
    *dt_cfl = damage_timestep_factor / dD_dt;
  }
}

/**
 * @brief Set the (symmetric) stress tensor by combining the deviatoric with the
 * pressure and applying damage.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void damage_compute_stress_tensor(
    struct sym_matrix *stress_tensor, const struct sym_matrix damaged_deviatoric_stress_tensor, const float pressure, const float damage) {

  *stress_tensor = damaged_deviatoric_stress_tensor;

  if (pressure < 0.f) {
    stress_tensor->xx -= (1.f - damage) * pressure;
    stress_tensor->yy -= (1.f - damage) * pressure;
    stress_tensor->zz -= (1.f - damage) * pressure;
  } else {
    stress_tensor->xx -= pressure;
    stress_tensor->yy -= pressure;
    stress_tensor->zz -= pressure;
  }
}

__attribute__((always_inline)) INLINE static void damage_reset_predicted_values(
    struct part *restrict p, const struct xpart *restrict xp) {
    
  strength_set_damage(p, xp->strength_data.damage_full);
  damage_set_tensile_damage(p, xp->strength_data.tensile_damage_full);
  damage_set_shear_damage(p, xp->strength_data.shear_damage_full);
}

/**
 * @brief Evolves particle damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void damage_evolve(
    struct part *restrict p, float *damage, float *tensile_damage, float *shear_damage,
    const struct sym_matrix stress_tensor, const struct sym_matrix deviatoric_stress_tensor, 
    const int mat_id, const float mass, const float density, const float u, const float yield_stress, const float dt_therm) {

    // ### note that time derivatives get calculated each time this gets called i.e. in all of kick-drift-kick
    // ### results are sensitive to how often these get recalculated so might need to do it each time like this
    // ### might be computationally expensive with recalculation of eigenvalues

  damage_tensile_evolve(p, tensile_damage, stress_tensor, mat_id, mass, density, *damage, dt_therm);

  damage_shear_evolve(p, shear_damage, deviatoric_stress_tensor, mat_id, density, u, yield_stress, dt_therm);

  *damage = fminf(*tensile_damage + *shear_damage, 1.f);
}

/**
 * @brief Evolves particle damage in the drift
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void damage_predict_evolve(
    struct part *restrict p, const struct sym_matrix stress_tensor, const struct sym_matrix deviatoric_stress_tensor, 
    const int mat_id, const float mass, const float density, const float u, const float yield_stress, const float dt_therm) {

  // Get damage parameters before they get evolved in time
  float damage = strength_get_damage(p);
  float tensile_damage = damage_get_tensile_damage(p);
  float shear_damage = damage_get_shear_damage(p);

  // Evolve damage parameters
  damage_evolve(p, &damage, &tensile_damage, &shear_damage, 
                  stress_tensor, deviatoric_stress_tensor, mat_id, mass, density, u, yield_stress, dt_therm);

  // Update damage particle properties
  strength_set_damage(p, damage);
  damage_set_tensile_damage(p, tensile_damage);
  damage_set_shear_damage(p, shear_damage);
}

/**
 * @brief Evolves particle damage in the kick
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void damage_kick_evolve(
    struct part *restrict p, struct xpart *restrict xp, const struct sym_matrix stress_tensor, const struct sym_matrix deviatoric_stress_tensor, 
    const int mat_id, const float mass, const float density, const float u, const float yield_stress, const float dt_therm) {

// Get damage parameters before they get evolved in time
  float damage = strength_get_damage_full(xp);
  float tensile_damage = damage_get_tensile_damage_full(xp);
  float shear_damage = damage_get_shear_damage_full(xp);

  // Evolve damage parameters
  damage_evolve(p, &damage, &tensile_damage, &shear_damage, 
                  stress_tensor, deviatoric_stress_tensor, mat_id, mass, density, u, yield_stress, dt_therm);

  // Update damage particle properties
  strength_set_damage_full(xp, damage);
  damage_set_tensile_damage_full(xp, tensile_damage);
  damage_set_shear_damage_full(xp, shear_damage);
}
    
/**
 * @brief Calculate time derivative of damage.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void damage_compute_dD_dt(
    struct part *restrict p, const struct sym_matrix stress_tensor, const struct sym_matrix deviatoric_stress_tensor, 
    const int mat_id, const float mass, const float density, const float u, const float yield_stress) {

  const float damage = strength_get_damage(p);
  const float tensile_damage = damage_get_tensile_damage(p);
  const float shear_damage = damage_get_shear_damage(p);

  float tensile_dD_dt = 0.f;
  damage_tensile_compute_dD_dt(p, &tensile_dD_dt, stress_tensor, mat_id, mass, density, damage, tensile_damage);
    
  float shear_dD_dt = 0.f;
  damage_shear_compute_dD_dt(&shear_dD_dt, deviatoric_stress_tensor, mat_id, density, u, yield_stress, shear_damage);

  p->strength_data.dD_dt = tensile_dD_dt + shear_dD_dt;
}

__attribute__((always_inline)) INLINE static void damage_first_init_part(
    struct part *restrict p, struct xpart *restrict xp) {

  strength_set_damage(p, 0.f);
  damage_set_tensile_damage(p, 0.f);
  damage_set_shear_damage(p, 0.f);
    
  strength_set_damage_full(xp, 0.f);
  damage_set_tensile_damage_full(xp, 0.f);
  damage_set_shear_damage_full(xp, 0.f);
}

#endif /* SWIFT_DAMAGE_H */
