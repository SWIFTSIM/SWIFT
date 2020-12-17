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
#ifndef LOGGER_LOGGER_PARTICLE_H
#define LOGGER_LOGGER_PARTICLE_H

/* Include the tools. */
#include "logger_tools.h"

/* Include the other local files. */
#include "logger_gravity.h"
#include "logger_header.h"
#include "logger_hydro.h"
#include "logger_stars.h"
#include "logger_time.h"

#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Defines the type of interpolation
 */
enum logger_reader_type {
  logger_reader_const, /* Constant interpolation. */
  logger_reader_lin,   /* Linear interpolation. */
};

size_t logger_particle_read_field(const struct logger_reader *reader,
                                  size_t offset, void *output,
                                  const int *local_to_global,
                                  const int local_count, int field,
                                  size_t *mask, size_t *h_offset);

enum logger_special_flags logger_particle_read_special_flag(
    const struct logger_reader *reader, size_t offset, size_t *mask,
    size_t *h_offset, int *data);

/**
 * @brief Generate the data for the special flags.
 *
 * @param flag The special flag to use.
 * @param data The data to write in the .
 */
INLINE static enum logger_special_flags logger_unpack_flags_and_data(
    uint32_t flag, int *data) {
  *data = flag & 0xFFFFFF;
  return flag >> (3 * 8);
}

#endif  // LOGGER_LOGGER_PARTICLE_H
