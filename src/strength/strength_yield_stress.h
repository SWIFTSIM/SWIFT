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
#ifndef SWIFT_STRENGTH_YIELD_STRESS_H
#define SWIFT_STRENGTH_YIELD_STRESS_H

/**
 * @file strength/strength_yield_stress.h
 * @brief Selects yield stress model based on configuration options.
 */

#if defined(STRENGTH_YIELD_STRESS_BENZ_ASPHAUG)
#include "yield_stress/yield_stress_ba95.h"
#elif defined(STRENGTH_YIELD_STRESS_COLLINS)
#include "yield_stress/yield_stress_collins04.h"
#else
#include "yield_stress/yield_stress_none.h"
#endif

#endif /* SWIFT_STRENGTH_YIELD_STRESS_H */
