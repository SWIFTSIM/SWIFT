/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_PRESSURE_FLOOR_IACT_H
#define SWIFT_PRESSURE_FLOOR_IACT_H

/**
 * @file src/pressure_floor_iact.h
 * @brief Branches between the different pressure floor iact.
 */

/* Config parameters. */
#include <config.h>

/* Import the right pressure floor definition */
#if defined(PRESSURE_FLOOR_NONE)
#include "./pressure_floor/none/pressure_floor_iact.h"
#elif defined(PRESSURE_FLOOR_GEAR)
#include "./pressure_floor/GEAR/pressure_floor_iact.h"
#else
#error "Invalid choice of pressure floor"
#endif

#endif /* SWIFT_PRESSURE_FLOOR_IACT_H */
