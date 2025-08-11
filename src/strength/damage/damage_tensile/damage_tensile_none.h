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
#ifndef SWIFT_DAMAGE_TENSILE_NONE_H
#define SWIFT_DAMAGE_TENSILE_NONE_H

/**
 * @file strength/damage/damage_tensile/damage_tensile_none.h
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
    struct part *restrict p, float *tensile_cbrtD_dt, float *number_of_activated_flaws, 
    struct sym_matrix deviatoric_stress_tensor, const float damage, const float density, const float u) {}

/**
 * @brief Evolves particle tensile damage
 *
 * @param p The particle to act upon
 */
__attribute__((always_inline)) INLINE static void evolve_damage_tensile(
    struct part *restrict p, float *tensile_damage, const float tensile_cbrtD_dt, 
    const float number_of_activated_flaws, const float dt_therm) {}

#endif /* SWIFT_DAMAGE_TENSILE_NONE_H */