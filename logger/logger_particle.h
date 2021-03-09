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
#include "logger_chemistry.h"
#include "logger_gravity.h"
#include "logger_header.h"
#include "logger_hydro.h"
#include "logger_parameters.h"
#include "logger_star_formation.h"
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
                                  const struct field_information *all_fields,
                                  const int all_fields_count,
                                  const int global_index, size_t *mask,
                                  size_t *h_offset);

void logger_particle_interpolate_field(
    const double t_before, const struct logger_field *restrict before,
    const double t_after, const struct logger_field *restrict after,
    void *restrict output, const double t,
    const struct field_information *field, enum part_type type,
    const struct logger_parameters *params);

enum logger_special_flags logger_particle_read_special_flag(
    const struct logger_reader *reader, size_t offset, size_t *mask,
    size_t *h_offset, int *data);

/**
 * @brief Extract the flag and its related data from the data inside a record.
 *
 * Flags and data are stored as a 32-bit unsigned int in which the lower 24 bits
 * contain the data, and the top 8 bits contain the flag value.

 * @param logfile_data The raw data taken from the logfile.
 * @param data (output) The data related to the flag.
 *
 * @return The flag extracted from the raw data.
 */

INLINE static enum logger_special_flags logger_unpack_flags_and_data(
    uint32_t logfile_data, int *flag_data) {
  const int num_data_bits = 24;
  *flag_data = logfile_data & ((1U << num_data_bits) - 1);
  return (enum logger_special_flags)(logfile_data >> num_data_bits);
}

#endif  // LOGGER_LOGGER_PARTICLE_H
