/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025   Jacob Kegerreis (j.kegerreis@imperial.ac.uk).
 *               2025   Thomas Sandnes (thomas.d.sandnes@durham.ac.uk).
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
#ifndef SWIFT_STRENGTH_H
#define SWIFT_STRENGTH_H

/* Config parameters. */
#include <config.h>

/* Import the right functions */
#ifdef MATERIAL_STRENGTH
#include "./hydro/REMIX/hydro_strength.h"
#include "./strength/strength_utilities.h"
#include "./strength/strength_utilities.h"
#include "./strength/strength_damage.h"
#include "./strength/strength_stress.h"
#include "./strength/strength_yield.h"
#endif /* MATERIAL_STRENGTH */

#endif /* SWIFT_STRENGTH_H */
