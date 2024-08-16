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
#ifndef SWIFT_PLANETARY_DAMAGE_H
#define SWIFT_PLANETARY_DAMAGE_H
#ifdef MATERIAL_STRENGTH

/**
 * @file Planetary/strength_damage.h
 * @brief REMIX implementation of SPH with material strength
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_kernels.h"
#include "hydro_parameters.h"
#include "hydro_strength.h"
#include "math.h"
#include "strength_utilities.h"


/**
 * @brief Set the (symmetric) stress tensor by combining the deviatoric with the
 * pressure and applying damage.
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static struct sym_matrix stress_tensor_damaged(
    struct sym_matrix deviatoric_stress_tensor, const float pressure, const float damage) {
    

    struct sym_matrix stress_tensor = deviatoric_stress_tensor;

    #if defined(STRENGTH_YIELD_BENZ_ASPHAUG) && defined(STRENGTH_DAMAGE)
      for (int i = 0; i < 6; i++) stress_tensor.elements[i] *= (1.f - damage);
    #endif /* STRENGTH_YIELD_BENZ_ASPHAUG && STRENGTH_DAMAGE */

    if (pressure < 0.f) {
      stress_tensor.xx -= (1.f - damage) * pressure;
      stress_tensor.yy -= (1.f - damage) * pressure;
      stress_tensor.zz -= (1.f - damage) * pressure;
    } else {
      stress_tensor.xx -= pressure;
      stress_tensor.yy -= pressure;
      stress_tensor.zz -= pressure;
    }

    return stress_tensor;
}

/**
 * @brief Adjust the yield stress depending on the damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float
adjust_yield_stress_by_damage(struct part *restrict p,
  const float yield_stress_intact, const float yield_stress_fully_damaged, const float damage) {

  float yield_stress = yield_stress_intact;

  #if defined(STRENGTH_DAMAGE)
        yield_stress = (1.f - damage) * yield_stress_intact +
                   damage * yield_stress_fully_damaged;
  #endif /* STRENGTH_DAMAGE */

  return yield_stress;
}


/**
 * @brief Calculates the rate of cbrt(damage) accumulated due to tension
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void calculate_tensile_cbrtD_dt(
    struct part *restrict p, float *tensile_cbrtD_dt, float *number_of_activated_flaws, 
    struct sym_matrix deviatoric_stress_tensor, const float damage, const float density, const float u) {

  #if defined(STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG)
  *tensile_cbrtD_dt = 0.f;  
  // Damage can only accumulate if in tension (this will have to change if we
  // add fracture under compression)

  // If there are no flaws, we can not accumulate damage
  if (p->number_of_flaws > 0) {
     
    // Compute current stress tensor to get max eigenvalue
    const float pressure =
       gas_pressure_from_internal_energy(density, u, p->mat_id);    
       
    struct sym_matrix stress_tensor = stress_tensor_damaged(deviatoric_stress_tensor, pressure, damage);
    float principal_stress_eigen[3];
    sym_matrix_compute_eigenvalues(principal_stress_eigen, stress_tensor); 

    // Find max eignenvalue of stress_tensor_sigma:
    float max_principal_stress = principal_stress_eigen[0];
    if (principal_stress_eigen[1] > max_principal_stress)
      max_principal_stress = principal_stress_eigen[1];
    if (principal_stress_eigen[2] > max_principal_stress)
      max_principal_stress = principal_stress_eigen[2];
     
    if (max_principal_stress > 0.f) {
      const float shear_mod = material_shear_mod(p->mat_id);
      const float bulk_mod = material_bulk_mod(p->mat_id);

      // tensile damage
      const float E = 9.f * bulk_mod * shear_mod / (3.f * bulk_mod + shear_mod);

      const float local_scalar_strain =
          max_principal_stress / ((1.f - damage) * E);

      *number_of_activated_flaws = 0;
      for (int i = 0; i < p->number_of_flaws; i++) {
        if (local_scalar_strain > p->activation_thresholds[i]) {
          *number_of_activated_flaws += 1;
        }
      }

      // Schafer
      const float longitudinal_wave_speed = sqrtf(
          (bulk_mod + (4.f / 3.f) * (1.f - damage) * shear_mod) / density);
      const float crack_velocity = 0.4f * longitudinal_wave_speed;

      *tensile_cbrtD_dt =
          *number_of_activated_flaws * crack_velocity * cbrtf(density / p->mass);
    }
  }
  #endif /* STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG */
}

/**
 * @brief Evolves particle tensile damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void evolve_damage_tensile(
    struct part *restrict p, float *tensile_damage, const float tensile_cbrtD_dt, 
    const float number_of_activated_flaws, const float dt_therm) {

  #if defined(STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG)
    float Delta_cbrtD = tensile_cbrtD_dt * dt_therm;

    const float max_cbrtD =
        cbrtf(number_of_activated_flaws / (float)p->number_of_flaws);
    const float max_Delta_cbrtD =
        max(0.f, max_cbrtD - cbrtf(*tensile_damage));

    if (Delta_cbrtD > max_Delta_cbrtD) Delta_cbrtD = max_Delta_cbrtD;

    const float evolved_D_cbrt = cbrtf(*tensile_damage) + Delta_cbrtD;
    *tensile_damage = powf(evolved_D_cbrt, 3.f);
  #endif /* STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG */
}

/**
 * @brief Calculates the rate of damage accumulated due to tension
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void calculate_shear_dD_dt(
    struct part *restrict p, float *shear_dD_dt, struct sym_matrix deviatoric_stress_tensor, const float yield_stress, const float density, const float u) {
   
  #if defined(STRENGTH_DAMAGE_SHEAR_COLLINS)
    *shear_dD_dt = 0.f;

    // ### Made some changes here. Compare with old version when testing
    
    const float J_2 = J_2_from_stress_tensor(deviatoric_stress_tensor);
    
    // See Collins for this.
    if (sqrtf(J_2) > yield_stress) {
      const float pressure =
          gas_pressure_from_internal_energy(density, u, p->mat_id);     
      // do I need this or can it happen when p is negative as well?
      if (pressure > 0) {
        // shear damage
        float brittle_to_ductile_pressure =
            material_brittle_to_ductile_pressure(p->mat_id);
        float brittle_to_plastic_pressure =
            material_brittle_to_plastic_pressure(p->mat_id);

        float plastic_strain_at_failure;
        if (pressure < brittle_to_ductile_pressure) {

          // Is this meant to be inear? maybe not clear in paper
          plastic_strain_at_failure =
              0.04f * (pressure / brittle_to_ductile_pressure) + 0.01f;

        } else if (pressure < brittle_to_plastic_pressure) {

          float slope = (0.1f - 0.05f) / (brittle_to_plastic_pressure -
                                        brittle_to_ductile_pressure);
          float intercept = brittle_to_ductile_pressure - slope * 0.05f;

          plastic_strain_at_failure =
              slope * pressure + intercept;

        } else {
          plastic_strain_at_failure = 1.f;
        }

        strain_rate_invariant = sqrtf(J_2);

        *shear_dD_dt = (strain_rate_invariant / plastic_strain_at_failure);
      }
    }
  #endif /* STRENGTH_DAMAGE_SHEAR_COLLINS */

}

/**
 * @brief Evolves particle shear damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void evolve_damage_shear(
    struct part *restrict p, float *shear_damage, const float  shear_dD_dt, const float dt_therm) {

  #if defined(STRENGTH_DAMAGE_SHEAR_COLLINS)
 
      *shear_damage += min(shear_dD_dt * dt_therm, 1.f);
    
  #endif /* STRENGTH_DAMAGE_SHEAR_COLLINS */
}

/**
 * @brief Evolves particle damage
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

#endif /* MATERIAL_STRENGTH */
#endif /* SWIFT_PLANETARY_DAMAGE_H */
