/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_IO_M1CLOSURE_H
#define SWIFT_RT_IO_M1CLOSURE_H

#include "io_properties.h"

/**
 * @file src/rt/M1closure/rt_io.h
 * @brief Main header file for no radiative transfer scheme IO routines.
 */

/**
 * @brief Creates additional output fields for the radiative transfer data of
 * hydro particles. Intended to be called from inside hydro_write_particles().
 */
INLINE static void rt_write_xparts(const struct xpart* xparts,
                                   struct io_props* list, int* num_fields) {}

/**
 * @brief Creates additional output fields for the radiative transfer data of
 * star particles. Intended to be called from inside hydro_write_particles().
 */
INLINE static void rt_write_stars(const struct spart* sparts,
                                  struct io_props* list, int* num_fields) {}
#endif /* SWIFT_RT_IO_M1CLOSURE_H */
