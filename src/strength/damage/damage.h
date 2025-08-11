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

__attribute__((always_inline)) INLINE static float strength_get_damage(struct part *restrict p) {

  return p->strength_data.damage;
}

__attribute__((always_inline)) INLINE static float strength_get_damage_full(const struct xpart *restrict xp) {

  return xp->strength_data.damage_full;
}

/**
 * @brief Set the (symmetric) stress tensor by combining the deviatoric with the
 * pressure and applying damage.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void adjust_stress_tensor_by_damage(struct part *restrict p,
    struct sym_matrix stress_tensor, struct sym_matrix deviatoric_stress_tensor, const float pressure, const float damage) {

  stress_tensor = yield_model_adjust_deviatoric_stress_tensor_by_damage(deviatoric_stress_tensor, damage);

  if (pressure < 0.f) {
    stress_tensor.xx -= (1.f - damage) * pressure;
    stress_tensor.yy -= (1.f - damage) * pressure;
    stress_tensor.zz -= (1.f - damage) * pressure;
  } else {
    stress_tensor.xx -= pressure;
    stress_tensor.yy -= pressure;
    stress_tensor.zz -= pressure;
  }
}

#include "damage_models.h"

__attribute__((always_inline)) INLINE static void damage_timestep(
    const struct part *restrict p, float dt_cfl) {

    const float damage_timestep_factor = 0.01f; // ### Hardcoded for now. Treat this similarly to CFL
    
  if (p->strength_data.dD_dt * dt_cfl > damage_timestep_factor) {
    dt_cfl = damage_timestep_factor / p->strength_data.dD_dt;
  }
}

__attribute__((always_inline)) INLINE static void hydro_reset_predicted_values_extra_damage(
    struct part *restrict p, const struct xpart *restrict xp) {
    
  p->strength_data.damage = xp->strength_data.damage_full;
  p->strength_data.tensile_damage = xp->strength_data.tensile_damage_full;
  p->strength_data.shear_damage = xp->strength_data.shear_damage_full;
}

/**
 * @brief Evolves particle damage in the drift
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void evolve_damage(
    struct part *restrict p, float *tensile_damage, float *shear_damage,
    float *damage, struct sym_matrix deviatoric_stress_tensor, const float yield_stress, 
    const float density, const float u, const float dt_therm) {

    // ### note that time derivatives get calculated each time this gets called i.e. in all of kick-drift-kick
    // ### results are sensitive to how often these get recalculated so might need to do it each time like this
    // ### might be computationally expensive with recalculation of eigenvalues

    float tensile_cbrtD_dt, number_of_activated_flaws;
    calculate_tensile_cbrtD_dt(p, &tensile_cbrtD_dt, &number_of_activated_flaws, 
                                        deviatoric_stress_tensor, *damage, density, u);
    evolve_damage_tensile(p, tensile_damage, tensile_cbrtD_dt, number_of_activated_flaws, dt_therm);


    float shear_dD_dt;
    calculate_shear_dD_dt(p, &shear_dD_dt, deviatoric_stress_tensor, yield_stress, density, u);
    evolve_damage_shear(p, shear_damage, shear_dD_dt, dt_therm);

    *damage = min(*tensile_damage + *shear_damage, 1.f);
}
    
/**
 * @brief Evolves particle damage in the drift
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_predict_evolve_damage(
    struct part *restrict p, const float yield_stress, const float density, const float u, const float dt_therm) {

    // ### note that time derivatives get calculated each time this gets called i.e. in all of kick-drift-kick
    // ### results are sensitive to how often these get recalculated so might need to do it each time like this
    // ### might be computationally expensive with recalculation of eigenvalues

   evolve_damage(p, &p->strength_data.tensile_damage, &p->strength_data.shear_damage, &p->strength_data.damage, p->strength_data.deviatoric_stress_tensor, yield_stress, density, u, dt_therm);
}

/**
 * @brief Evolves particle damage in the kick
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void hydro_kick_evolve_damage(
    struct part *restrict p, struct xpart *restrict xp, const float yield_stress, const float density, const float u, const float dt_therm) {

    // ### note that time derivatives get calculated each time this gets called i.e. in all of kick-drift-kick
    // ### results are sensitive to how often these get recalculated so might need to do it each time like this
    // ### might be computationally expensive with recalculation of eigenvalues

    evolve_damage(p, &xp->strength_data.tensile_damage_full, &xp->strength_data.shear_damage_full, &xp->strength_data.damage_full, xp->strength_data.deviatoric_stress_tensor_full, yield_stress, density, u,  dt_therm);
}

/**
 * @brief Calculate time derivative of damage.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void update_dD_dt(
    struct part *restrict p, const int phase_state, const float density, const float u) {

  const float damage = p->strength_data.damage;
  const float yield_stress = compute_yield_stress_damage(p, phase_state, density, u, damage);

  p->strength_data.dD_dt = 0.f;
    
  float tensile_cbrtD_dt = 0.f;
  float number_of_activated_flaws = 0;
  calculate_tensile_cbrtD_dt(p, &tensile_cbrtD_dt, &number_of_activated_flaws, 
                                        p->strength_data.deviatoric_stress_tensor, damage, density, u);
  if (p->strength_data.tensile_damage < number_of_activated_flaws / (float)p->strength_data.number_of_flaws) {
    // Chain rule d(D^(1/3))/dt  = d(D^(1/3))/dD * dD/dt
    p->strength_data.dD_dt += 3.f * powf(p->strength_data.tensile_damage, 2.f / 3.f) * tensile_cbrtD_dt;
  }
    
  float shear_dD_dt = 0.f;
  calculate_shear_dD_dt(p, &shear_dD_dt, p->strength_data.deviatoric_stress_tensor, yield_stress, density, u);

  if (p->strength_data.shear_damage < 1.f) {
    p->strength_data.dD_dt += shear_dD_dt;
  } 
}

__attribute__((always_inline)) INLINE static void hydro_first_init_part_damage(
    struct part *restrict p, struct xpart *restrict xp) {

  p->strength_data.damage = 0.f;
  p->strength_data.tensile_damage = 0.f;
  p->strength_data.shear_damage = 0.f;

  xp->strength_data.damage_full = 0.f;
  xp->strength_data.tensile_damage_full = 0.f;
  xp->strength_data.shear_damage_full = 0.f;
}

#endif /* SWIFT_DAMAGE_H */
