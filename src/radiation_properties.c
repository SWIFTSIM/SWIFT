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

/* This object's header. */
#include "radiation_properties.h"

/* Standard headers */
#include <float.h>
#include <math.h>

/* Local headers. */
#include "common_io.h"
#include "error.h"

void radiation_props_init(struct rad_props *p, struct swift_params *params) {}

void radiation_props_print(const struct rad_props *p) {}

#if defined(HAVE_HDF5)
void radiation_props_print_snapshot(hid_t h_grpgrav,
                                    const struct rad_props *p) {}
#endif

/**
 * @brief Write a rad_props struct to the given FILE as a stream of bytes.
 *
 * @param p the struct
 * @param stream the file stream
 */
void radiation_props_struct_dump(const struct rad_props *p, FILE *stream) {
  restart_write_blocks((void *)p, sizeof(struct rad_props), 1, stream,
                       "radiation", "rad props");
}

/**
 * @brief Restore a rad_props struct from the given FILE as a stream of
 * bytes.
 *
 * @param p the struct
 * @param stream the file stream
 */
void radiation_props_struct_restore(struct rad_props *p, FILE *stream) {
  restart_read_blocks((void *)p, sizeof(struct rad_props), 1, stream, NULL,
                      "rad props");
}
