/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2022 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_STAR_FORMATION_DEBUG_H
#define SWIFT_STAR_FORMATION_DEBUG_H

/* Config parameters. */
#include "../config.h"

/* Import the debug routines of the right star formation definition */
#if defined(STAR_FORMATION_NONE)
#include "./star_formation/none/star_formation_debug.h"
#elif defined(STAR_FORMATION_QLA)
#include "./star_formation/QLA/star_formation_debug.h"
#elif defined(STAR_FORMATION_EAGLE)
#include "./star_formation/EAGLE/star_formation_debug.h"
#elif defined(STAR_FORMATION_GEAR)
#include "./star_formation/GEAR/star_formation_debug.h"
#else
#error "Invalid choice of star formation model."
#endif

#endif /* SWIFT_STAR_FORMATION_DEBUG_H */
