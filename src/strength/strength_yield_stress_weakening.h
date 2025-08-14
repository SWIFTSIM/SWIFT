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
#ifndef SWIFT_STRENGTH_YIELD_STRESS_WEAKENING_H
#define SWIFT_STRENGTH_YIELD_STRESS_WEAKENING_H

/**
 * @file strength/strength_yield_stress_weakening.h
 * Methods for applying weakening to the yield stress.
 */

#include "equation_of_state.h"
#include "math.h"

/**
 * @brief Apply thermal weakening to the yield stress.
 *
 * Weaken the material's yield stress as its temperature approaches the
 * melting temperature. This method is presented in Ohnaka1995, and
 * is used in, e.g., Collins+2004 and Emsenhuber+2018. Note that we do not
 * set the yield stress to 0 here when the temperature exceeds the melting
 * temperature; this is handled elsewhere based on the material's phase state.
 *
 * @param Y the yield stress to be weakened.
 * @param mat_id the material ID.
 * @param density the density.
 * @param u the specific internal energy.
 */
__attribute__((always_inline)) INLINE static void
yield_weakening_apply_temperature_to_yield_stress(float *Y, const int mat_id,
                                   const float density, const float u) {
#ifdef STRENGTH_YIELD_STRESS_WEAKENING_THERMAL
  const float temperature =
        gas_temperature_from_internal_energy(density, u, mat_id);
  const float xi = method_yield_thermal_soft_xi();
  const float T_melt = material_T_melt(mat_id);
    
  *Y *= tanhf(xi * (T_melt / temperature - 1.f));
#endif /* STRENGTH_YIELD_STRESS_WEAKENING_THERMAL */
}

/**
 * @brief Apply density weakening to the yield stress.
 *
 * Weaken the material's yield stress as its density decreases below a certain
 * point. This method is discussed in Luther+2022, where it is described in the
 * context of its implementation in iSALE.
 *
 * @param Y the yield stress to be weakened.
 * @param mat_id the material ID.
 * @param density the density.
 */
__attribute__((always_inline)) INLINE static void
yield_weakening_apply_density_to_yield_stress(float *Y, const int mat_id,
                               const float density) {
#ifdef STRENGTH_YIELD_STRESS_WEAKENING_DENSITY
  const float rho_0 = material_rho_0(mat_id);
  const float a = method_yield_density_soft_mult_param();
  const float rho_weak = a * rho_0;
    
  if (density < rho_weak) {
    const float b = method_yield_density_soft_pow_param();
    *Y *= powf(density / rho_weak, b);
  }
#endif /* STRENGTH_YIELD_STRESS_WEAKENING_DENSITY */
}

#endif /* SWIFT_STRENGTH_YIELD_STRESS_WEAKENING_H */