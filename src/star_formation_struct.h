/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_STAR_FORMATION_STRUCT_H
#define SWIFT_STAR_FORMATION_STRUCT_H

/**
 * @file src/star_formation_struct.h
 * @brief Branches between the different particle data SF tracers
 */

/* Config parameters. */
#include "../config.h"

/* Import the right cooling definition */
#if defined(STAR_FORMATION_NONE)
#include "./star_formation/none/star_formation_struct.h"
#elif defined(STAR_FORMATION_EAGLE)
#include "./star_formation/EAGLE/star_formation_struct.h"
#else
#error "Invalid choice of star formation structure."
#endif

#endif /* SWIFT_STAR_FORMATION_STRUCT_H */
