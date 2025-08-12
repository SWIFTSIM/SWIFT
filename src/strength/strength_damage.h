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
#ifndef SWIFT_STRENGTH_DAMAGE_H
#define SWIFT_STRENGTH_DAMAGE_H

/**
 * @file strength/strength_damage.h
 * @brief Selects whether using damage models or not on configuration options.
 */

#if defined(STRENGTH_DAMAGE)

#if defined(STRENGTH_DAMAGE_SHEAR_COLLINS)
#include "damage/damage_shear/damage_shear_collins04.h"
#else
#include "damage/damage_shear/damage_shear_none.h"
#endif /* STRENGTH_DAMAGE_SHEAR_COLLINS */

#if defined(STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG)
#include "damage/damage_tensile/damage_tensile_ba94.h"
#else
#include "damage/damage_tensile/damage_tensile_none.h"
#endif /* STRENGTH_DAMAGE_TENSILE_BENZ_ASPHAUG */

#include "damage/damage_core.h"
#else
#include "damage/damage_none.h"
#endif

#endif /* SWIFT_STRENGTH_DAMAGE_H */
