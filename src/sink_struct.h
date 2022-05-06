/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Yves Revaz (yves.revaz@epfl.ch)
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
#ifndef SWIFT_SINK_STRUCT_H
#define SWIFT_SINK_STRUCT_H

/**
 * @file src/sink_struct.h
 * @brief Branches between the different sink functions.
 */

/* Config parameters. */
#include "../config.h"
#include "inline.h"

/* Import the right black holes definition */
#if defined(SINK_NONE)
#include "./sink/Default/sink_struct.h"
#elif defined(SINK_GEAR)
#include "./sink/GEAR/sink_struct.h"
#else
#error "Invalid choice of sink model."
#endif

#endif /* SWIFT_SINK_STRUCT_H */
