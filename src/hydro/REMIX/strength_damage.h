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
 * @brief Compute the effective pressure from the damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float
effective_pressure_from_damage(struct part *restrict p, const float pressure) {

  float effective_pressure = pressure;

  #ifdef STRENGTH_DAMAGE
    // See schafer 2016 for this...
    if (pressure < 0.f) {
      effective_pressure *= (1.f - p->damage);
    }
  #endif /* STRENGTH_DAMAGE */

  return effective_pressure;
}

/**
 * @brief Adjust the yield stress depending on the damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float
adjust_yield_stress_by_damage(struct part *restrict p,
  const float yield_stress_intact, const float yield_stress_damaged) {

  float yield_stress = yield_stress_intact;

  #if defined(STRENGTH_YIELD_COLLINS) && defined(STRENGTH_DAMAGE)
        yield_stress = (1.f - p->damage) * yield_stress_intact +
                   p->damage * yield_stress_damaged;
  #endif /* STRENGTH_YIELD_COLLINS && STRENGTH_DAMAGE */

  return yield_stress;
}

/**
 * @brief Evolves particle tensile damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void evolve_damage_tensile(
    struct part *restrict p, float *tensile_damage, const float dt_therm) {

  #if defined(STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG)
  // Find max eignenvalue of stress_tensor_sigma:
  float max_principal_stress = p->principal_stress_eigen[0];
  if (p->principal_stress_eigen[1] > max_principal_stress)
    max_principal_stress = p->principal_stress_eigen[1];
  if (p->principal_stress_eigen[2] > max_principal_stress)
    max_principal_stress = p->principal_stress_eigen[2];

  float cbrtD_dt = 0.f;
  // Damage can only accumulate if in tension (this will have to change if we
  // add fracture under compression)

  // If there are no flaws, we can not accumulate damage
  if (p->number_of_flaws > 0) {
    if (max_principal_stress > 0.f) {
      const float shear_mod = material_shear_mod(p->mat_id);
      const float bulk_mod = material_bulk_mod(p->mat_id);

      // tensile damage
      const float E = 9.f * bulk_mod * shear_mod / (3.f * bulk_mod + shear_mod);

      const float local_scalar_strain =
          max_principal_stress / ((1.f - p->damage) * E);

      int number_of_activated_flaws = 0;
      for (int i = 0; i < p->number_of_flaws; i++) {
        if (local_scalar_strain > p->activation_thresholds[i]) {
          number_of_activated_flaws += 1;
        }
      }

      // Schafer
      const float longitudinal_wave_speed = sqrtf(
          (bulk_mod + (4.f / 3.f) * (1.f - p->damage) * shear_mod) / p->rho);
      const float crack_velocity = 0.4f * longitudinal_wave_speed;

      cbrtD_dt =
          number_of_activated_flaws * crack_velocity * cbrtf(p->rho / p->mass);

      float Delta_cbrtD = cbrtD_dt * dt_therm;

      const float max_cbrtD =
          cbrtf(number_of_activated_flaws / (float)p->number_of_flaws);
      const float max_Delta_cbrtD =
          max(0.f, max_cbrtD - cbrtf(*tensile_damage));

      if (Delta_cbrtD > max_Delta_cbrtD) Delta_cbrtD = max_Delta_cbrtD;

      const float evolved_D_cbrt = cbrtf(*tensile_damage) + Delta_cbrtD;
      *tensile_damage = powf(evolved_D_cbrt, 3.f);
    }
  }
  #endif /* STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG */
}

/**
 * @brief Evolves particle shear damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void evolve_damage_shear(
    struct part *restrict p, float *shear_damage, const float pressure,
    const float dt_therm) {

  #if defined(STRENGTH_DAMAGE_SHEAR_COLLINS)
  float J_2 = J_2_from_stress_tensor(&p->deviatoric_stress_tensor);

  // See Collins for this.
  if (sqrtf(J_2) > p->yield_stress) {
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

      // ### The next part feels like a recalculation of J_2. Double check this
      // Do I need to have rotation terms here?
      float strain_rate_tensor[3][3];
      calculate_strain_rate_tensor(p, strain_rate_tensor);

      float deviatoric_strain_rate_tensor[3][3];
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          deviatoric_strain_rate_tensor[i][j] = strain_rate_tensor[i][j];
          if (i == j) {
            deviatoric_strain_rate_tensor[i][j] -=
                (strain_rate_tensor[0][0] + strain_rate_tensor[1][1] +
                 strain_rate_tensor[2][2]) /
                3.f;
          }
        }
      }

      float strain_rate_invariant = 0.f;
      for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
          strain_rate_invariant += 0.5f * deviatoric_strain_rate_tensor[i][j] *
                                   deviatoric_strain_rate_tensor[i][j];
        }
      }

      strain_rate_invariant = sqrtf(strain_rate_invariant);

      // can this be negative?
      *shear_damage +=
          max(0.f, (strain_rate_invariant /
                    plastic_strain_at_failure) *
                       dt_therm);

      if (*shear_damage > 1.f) *shear_damage = 1.f;
    }
  }
  #endif /* STRENGTH_DAMAGE_SHEAR_COLLINS */
}

/**
 * @brief Evolves particle damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void evolve_damage(
    struct part *restrict p, const float pressure, const float dt_therm) {

  #if defined(STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG)
    evolve_damage_tensile(p, &p->tensile_damage, dt_therm);
  #endif

  #if defined(STRENGTH_DAMAGE_SHEAR_COLLINS)
    evolve_damage_shear(p, &p->shear_damage, pressure, dt_therm);
  #endif

  #if defined(STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG) && defined(STRENGTH_DAMAGE_SHEAR_COLLINS)
    p->damage = min(1.f, p->tensile_damage + p->shear_damage);
  #elif defined(STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG)
    p->damage = p->tensile_damage;
  #elif defined(STRENGTH_DAMAGE_SHEAR_COLLINS)
    p->damage = p->shear_damage;
  #endif
}

#endif /* MATERIAL_STRENGTH */
#endif /* SWIFT_PLANETARY_DAMAGE_H */
