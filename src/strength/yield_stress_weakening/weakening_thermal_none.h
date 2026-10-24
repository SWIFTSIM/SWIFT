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
#ifndef SWIFT_WEAKENING_THERMAL_NONE_H
#define SWIFT_WEAKENING_THERMAL_NONE_H

/**
 * @file strength/yield_stress_weakening/weakening_thermal_none.h
 * @brief No thermal weakening of the yield stress.
 */

#include "math.h"

/**
 * @brief Apply thermal weakening to the yield stress.
 *
 * Empty function when configuring without thermal weakening.
 *
 * @param Y The yield stress to be weakened.
 * @param density The density.
 * @param u The specific internal energy.
 * @param mat_id The material ID.
 * @param mp The material parameters.
 */
__attribute__((always_inline)) INLINE static void
yield_weakening_apply_temperature_to_yield_stress(float *Y, const float density,
                                                  const float u, const int mat_id,
                                                  const struct mat_params *mp) {}

#endif /* SWIFT_WEAKENING_THERMAL_NONE_H */
