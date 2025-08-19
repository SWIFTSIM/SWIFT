/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_DYNAMICAL_FRICTION_STRUCT_H
#define SWIFT_DYNAMICAL_FRICTION_STRUCT_H

/**
 * @file src/dynamical_friction_struct.h
 * @brief Branches between the different dynamical friction functions.
 */

/* Config parameters. */
#include <config.h>

/* Import the right dynamical friction definition */
#if defined(DYNAMICAL_FRICTION_NONE)
#include "./dynamical_friction/none/dynamical_friction_struct.h"
#elif defined(DYNAMICAL_FRICTION_BASIC)
#include "./dynamical_friction/Basic/dynamical_friction_struct.h"
#else
#error "Invalid choice of dynamical friction function."
#endif

#endif /* SWIFT_DYNAMICAL_FRICTION_STRUCT_H */
