/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_TRACERS_IO_H
#define SWIFT_TRACERS_IO_H

/**
 * @file src/tracers_io.h
 * @brief Branches between the different particle data tracers
 */

/* Config parameters. */
#include "../config.h"

/* Import the right cooling definition */
#if defined(TRACERS_NONE)
#include "./tracers/none/tracers_io.h"
#elif defined(TRACERS_EAGLE)
#include "./tracers/EAGLE/tracers_io.h"
#else
#error "Invalid choice of tracers."
#endif

#endif /* SWIFT_TRACERS_IO_H */
