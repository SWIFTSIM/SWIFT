/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_RADIATION_PROPERTIES
#define SWIFT_RADIATION_PROPERTIES

/* Config parameters. */
#include "../config.h"

#if defined(HAVE_HDF5)
#include <hdf5.h>
#endif

/* Local includes. */
#include "restart.h"

/* Forward declarations */
struct swift_params;

/**
 * @brief Contains all the constants and parameters of the radiation transfer
 * scheme.
 */
struct rad_props {};

void radiation_props_print(const struct rad_props *p);
void radiation_props_init(struct rad_props *p, struct swift_params *params);

#if defined(HAVE_HDF5)
void radiation_props_print_snapshot(hid_t h_grpsph, const struct rad_props *p);
#endif

/* Dump/restore. */
void radiation_props_struct_dump(const struct rad_props *p, FILE *stream);
void radiation_props_struct_restore(struct rad_props *p, FILE *stream);

#endif /* SWIFT_RADIATION_PROPERTIES */
