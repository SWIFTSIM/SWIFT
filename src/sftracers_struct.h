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
#ifndef SWIFT_SFTRACERS_STRUCT_H
#define SWIFT_SFTRACERS_STRUCT_H

/**
 * @file src/tracers_struct.h
 * @brief Branches between the different particle data tracers
 */

/* Config parameters. */
#include "../config.h"

/* Import the right cooling definition */
#if defined(SFTRACERS_NONE)
#include "./sftracers/none/sftracers_struct.h"
#elif defined(SFTRACERS_EAGLE)
#include "./sftracers/EAGLE/sftracers_struct.h"
#else
#error "Invalid choice of star formation tracers."
#endif

#endif /* SWIFT_SFTRACERS_STRUCT_H */
