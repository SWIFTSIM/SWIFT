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
#ifndef SWIFT_FEEDBACK_STRUCT_H
#define SWIFT_FEEDBACK_STRUCT_H

/**
 * @file src/feedback_struct.h
 * @brief Branches between the different feedback functions.
 */

/* Config parameters. */
#include <config.h>

/* Import the right feedback definition */
#if defined(FEEDBACK_NONE)
#include "./feedback/none/feedback_struct.h"
#elif defined(FEEDBACK_EAGLE_THERMAL)
#include "./feedback/EAGLE_thermal/feedback_struct.h"
#elif defined(FEEDBACK_EAGLE_KINETIC)
#include "./feedback/EAGLE_kinetic/feedback_struct.h"
#elif defined(FEEDBACK_GEAR)
#include "./feedback/GEAR/feedback_struct.h"
#elif defined(FEEDBACK_GEAR_MECHANICAL)
#include "./feedback/GEAR_mechanical/feedback_struct.h"
#elif defined(FEEDBACK_AGORA)
#include "./feedback/AGORA/feedback_struct.h"
#else
#error "Invalid choice of feedback function."
#endif

#endif /* SWIFT_FEEDBACK_STRUCT_H */
