/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 John Helly (j.c.helly@durham.ac.uk)
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
#ifndef SWIFT_HBT_INTERFACE_H
#define SWIFT_HBT_INTERFACE_H

/* Config parameters. */
#include "../config.h"

/* Forward declaration */
struct engine;

#ifdef HAVE_HBT

/**
 * @brief Initialise HBT library
 */
void hbt_init(struct engine *e);

/**
 * @brief Invoke the HBT halo finder
 *
 * @param output_nr Index of the current output
 *
 */
void hbt_invoke(struct engine *e, const int output_nr);

/**
 * @brief Free resources used by HBT library
 */
void hbt_free(void);

#endif /* HAVE_HBT */

#endif /* SWIFT_HBT_INTERFACE_H */
