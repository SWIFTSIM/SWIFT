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
#ifndef SWIFT_STRENGTH_YIELD_STRESS_SOFTENING_H
#define SWIFT_STRENGTH_YIELD_STRESS_SOFTENING_H

/**
 * @file strength/strength_yield_stress_softening.h
 */

#include "equation_of_state.h"
#include "hydro_parameters.h"
#include "math.h"
#include "symmetric_matrix.h"

/**
 * @brief Adjust the yield stress depending on the density
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float
yield_softening_apply_temperature_to_yield_stress(const float Y, const int mat_id,
                                   const float density, const float u) {

#ifdef STRENGTH_YIELD_STRESS_SOFTENING_THERMAL
  const float temperature =
        gas_temperature_from_internal_energy(density, u, mat_id);
  //  This is from Emsenhuber+2017
  //  See Senft+Stewart2007 for why this comes here and not just for intact
  const float xi = method_yield_thermal_soft_xi();
  const float T_melt = material_T_melt(mat_id);
    
  return Y * tanhf(xi * (T_melt / temperature - 1.f));
#else
  return Y;
#endif /* STRENGTH_YIELD_STRESS_SOFTENING_THERMAL */


}

/**
 * @brief Adjust the yield stress depending on the density
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static float
yield_softening_apply_density_to_yield_stress(const float Y, const int mat_id,
                               const float density) {

#ifdef STRENGTH_YIELD_STRESS_SOFTENING_DENSITY
  const float rho_0 = material_rho_0(mat_id);
  const float a = method_yield_density_soft_mult_param();
  const float rho_weak = a * rho_0;
    
  if (density < rho_weak) {
    const float b = method_yield_density_soft_pow_param();
    return Y * powf(density / rho_weak, b);
  } else {
    return Y;
  }
#else
  return Y;
#endif /* STRENGTH_YIELD_STRESS_SOFTENING_DENSITY */
}

#endif /* SWIFT_STRENGTH_YIELD_STRESS_SOFTENING_H */