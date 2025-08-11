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
#ifndef SWIFT_DAMAGE_SHEAR_COLLINS04_H
#define SWIFT_DAMAGE_SHEAR_COLLINS04_H

/**
 * @file strength/damage/damage_shear/damage_shear_collins04.h
 */

#include "const.h"
#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"


/**
 * @brief Calculates the rate of damage accumulated due to tension
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void calculate_shear_dD_dt(
    struct part *restrict p, float *shear_dD_dt, struct sym_matrix deviatoric_stress_tensor, const float yield_stress, const float density, const float u) {
   
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
}

/**
 * @brief Evolves particle shear damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void evolve_damage_shear(
    struct part *restrict p, float *shear_damage, const float  shear_dD_dt, const float dt_therm) {
 
  *shear_damage += min(shear_dD_dt * dt_therm, 1.f);
}

#endif /* SWIFT_DAMAGE_SHEAR_COLLINS04_H */