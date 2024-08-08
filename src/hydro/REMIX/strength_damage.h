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
#include "math.h"


/**
 * @brief Evolves particle tensile damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void set_effective_pressure_by_damage(
    const float *effective_pressure, const float pressure, const float damage) {

#if defined(STRENGTH_DAMAGE_###)
    effective_pressure = pressure;

    // See schafer 2016 for this...
    if (pressure < 0.f) {
        effective_pressure *= (1.f - p->damage);
    }
#endif
}

/**
 * @brief Evolves particle tensile damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void evolve_damage_tensile(
    struct part *restrict p, const float *tensile_damage, const float dt_therm) {

#if defined(STRENGTH_DAMAGE_###)
    // Find max eignenvalue of stress_tensor_sigma:
    float eigenval1, eigenval2, eigenval3;
    compute_eigenvalues_symmetric_3x3(&eigenval1, &eigenval2, &eigenval3, p->stress_tensor_sigma);
    float max_stress_tensor_sigma = eigenval1;
    if (eigenval2 > max_stress_tensor_sigma)
        max_stress_tensor_sigma = eigenval2;
    if (eigenval3 > max_stress_tensor_sigma)
        max_stress_tensor_sigma = eigenval3;

    p->cbrtD_dt = 0.f;
    // Damage can only accumulate if in tension (this will have to change if we
    // add fracture under compression)

    // If there are no flaws, we can not accumulate damage
    if(p->number_of_flaws > 0) {
        if (max_stress_tensor_sigma > 0.f) {
            const float shear_mod = material_shear_mod(p->mat_id);
            const float bulk_mod = material_bulk_mod(p->mat_id);

            // tensile damage
            const float E = 9.f * bulk_mod * shear_mod / (3.f * bulk_mod + shear_mod);

            p->local_scalar_strain = max_stress_tensor_sigma / ((1.f - p->damage) * E);

            int number_of_activated_flaws = 0;
            for (int i = 0; i < p->number_of_flaws; i++) {
                if (p->local_scalar_strain > p->activation_thresholds_epsilon_act_ij[i]){
                    number_of_activated_flaws += 1;
                }
            }
            p->number_of_activated_flaws = max(
                p->number_of_activated_flaws, number_of_activated_flaws);

            // Schafer
            const float longitudinal_wave_speed =
                sqrtf((bulk_mod + (4.f / 3.f) * (1.f - p->damage) * shear_mod) / p->rho);
            const float crack_velocity = 0.4f * longitudinal_wave_speed;

            p->cbrtD_dt = number_of_activated_flaws * crack_velocity * cbrtf(p->rho / p->mass);

            float Delta_cbrtD = p->cbrtD_dt * dt_therm;

            const float max_cbrtD = cbrtf(number_of_activated_flaws / (float)p->number_of_flaws);
            const float max_Delta_cbrtD = max(0.f, max_cbrtD - cbrtf(tensile_damage));

            if (Delta_cbrtD > max_Delta_cbrtD)
                Delta_cbrtD = max_Delta_cbrtD;

            const float evolved_D_cbrt = cbrtf(tensile_damage) + Delta_cbrtD;
            tensile_damage = powf(evolved_D_cbrt, 3.f);
        }
    }
#endif
}


/**
 * @brief Evolves particle shear damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void evolve_damage_shear(
    struct part *restrict p, const float pressure, const float dt_therm) {

#if defined(STRENGTH_DAMAGE_###)
    float J_2;
    compute_stress_tensor_J_2(&J_2, deviatoric_stress_tensor);

    // See Collins for this.
    if (sqrtf(J_2) > p->yield_stress_Y) {
      // do I need this or can it happen when p is negative as well?
      if (pressure > 0) {
        // shear damage
        float plastic_strain_at_the_point_of_failure_epsilon_f;
        if (pressure < p->brittle_to_ductile_transition_pressure) {

          // Is this meant to be inear? maybe not clear in paper
          plastic_strain_at_the_point_of_failure_epsilon_f =
              0.04f * (pressure / p->brittle_to_ductile_transition_pressure) +
              0.01f;

        } else if (pressure < p->brittle_to_plastic_transition_pressure) {

          float slope =
              (0.1f - 0.05f) / (p->brittle_to_plastic_transition_pressure -
                                p->brittle_to_ductile_transition_pressure);
          float intercept =
              p->brittle_to_ductile_transition_pressure - slope * 0.05f;

          plastic_strain_at_the_point_of_failure_epsilon_f =
              slope * pressure + intercept;

        } else {
          plastic_strain_at_the_point_of_failure_epsilon_f = 1.f;
        }

        // Do I need to have rotation terms here?
        float deviatoric_strain_rate_tensor[3][3];
        for (int i = 0; i < 3; i++) {
          for (int j = 0; j < 3; j++) {
            deviatoric_strain_rate_tensor[i][j] =
                p->strain_rate_tensor_epsilon_dot[i][j];
            if (i == j) {
              deviatoric_strain_rate_tensor[i][j] -=
                  (p->strain_rate_tensor_epsilon_dot[0][0] +
                   p->strain_rate_tensor_epsilon_dot[1][1] +
                   p->strain_rate_tensor_epsilon_dot[2][2]) /
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
        shear_damage +=
            max(0.f, (strain_rate_invariant /
                      plastic_strain_at_the_point_of_failure_epsilon_f) *
                         dt_therm);

        if (shear_damage > 1.f)
            shear_damage = 1.f;
      }
    }
#endif
}

/**
 * @brief Evolves particle damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void evolve_damage(
    struct part *restrict p, const float pressure, const float dt_therm) {

  // ...
  float tensile_damage = p->tensile_damage;
  float shear_damage = p->shear_damage;

  // ...
  evolve_damage_tensile(p, &tensile_damage, dt_therm);
  evolve_damage_shear(p, &shear_damage, pressure, dt_therm);

  // total damage
  p->damage = min(1.f, tensile_damage + shear_damage);
}

#endif /* MATERIAL_STRENGTH */
#endif /* SWIFT_PLANETARY_DAMAGE_H */
