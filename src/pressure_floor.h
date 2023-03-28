/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_PRESSURE_FLOOR_H
#define SWIFT_PRESSURE_FLOOR_H

/**
 * @file src/pressure_floor.h
 * @brief Branches between the different pressure floor models
 */

/* Config parameters. */
#include <config.h>

/* Local includes */
#include "common_io.h"
#include "error.h"
#include "inline.h"

/* Check if pressure floor is implemented in hydro */
#ifndef PRESSURE_FLOOR_NONE
#if defined(GADGET2_SPH) || defined(HOPKINS_PU_SPH) || defined(SPHENIX_SPH) || \
    defined(GASOLINE_SPH)
/* Implemented */
#else
#error Pressure floor not implemented with this hydro scheme
#endif
#endif

/* Check if the pressure floor is implemented in the stellar feedback */
#ifdef PRESSURE_FLOOR_GEAR
#if defined(FEEDBACK_EAGLE_THERMAL) || defined(FEEDBACK_EAGLE_KINETIC)
#error Pressure floor not implemented in this stellar feedback scheme
#endif
#endif

/* Check if the pressure floor is implemented in the AGN feedback */
#ifdef PRESSURE_FLOOR_GEAR
#if defined(BLACK_HOLES_EAGLE) || defined(BLACK_HOLES_SPIN_JET)
#error Pressure floor not implemented in this AGN feedback scheme
#endif
#endif

/* Import the right pressure floor definition */
#if defined(PRESSURE_FLOOR_NONE)
#include "./pressure_floor/none/pressure_floor.h"
#elif defined(PRESSURE_FLOOR_GEAR)
#include "./pressure_floor/GEAR/pressure_floor.h"
#endif

#endif /* SWIFT_PRESSURE_FLOOR_H */
