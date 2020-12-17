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
#include "logger_particle.h"

#include "logger_header.h"
#include "logger_loader_io.h"
#include "logger_reader.h"
#include "logger_time.h"
#include "logger_tools.h"

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

/**
 * @brief Read a particle (of any type) record in the log file.
 *
 * @param reader The #logger_reader.
 * @param offset offset of the record to read.
 * @param output Buffer for the requested field..
 * @param local_to_global The list converting local to global id for the fields
 * (e.g. gravity_logger_local_to_global)
 * @param local_count The number of element in logger_mask_id (e.g.
 * gravity_logger_field_count)
 * @param field The request field (in local indexing).
 * @param mask (out) The mask of the record.
 * @param h_offset (out) Difference of offset with the next record.
 *
 * @return position after the record.
 */
__attribute__((always_inline)) INLINE size_t logger_particle_read_field(
    const struct logger_reader *reader, size_t offset, void *output,
    const int *local_to_global, const int local_count, int field, size_t *mask,
    size_t *h_offset) {

  /* Get a few pointers. */
  const struct header *h = &reader->log.header;
  void *map = reader->log.log.map;

  *mask = 0;
  *h_offset = 0;

  /* Read the record's mask. */
  map = logger_loader_io_read_mask(h, (char *)map + offset, mask, h_offset);

  /* Check that the mask is within the limits. */
  if (*mask > (unsigned int)(1 << h->masks_count)) {
    error_python("Found an unexpected mask %zi", *mask);
  }

  /* Check if it is not a time record. */
  if (*mask == h->timestamp_mask) {
    error_python("Cannot read a particle from timestep record.");
  }

  /* Skip the special field. */
  if (*mask & h->masks[logger_index_special_flags].mask) {
    map += h->masks[logger_index_special_flags].size;
  }

  /* Read the record and copy it to a particle. */
  for (int local = 0; local < local_count; local++) {
    const int global = local_to_global[local];

    /* Is the mask present? */
    if (!(*mask & h->masks[global].mask)) {
      continue;
    }

    if (field == local) {
      /* Read the data. */
      map = logger_loader_io_read_data(map, h->masks[global].size, output);
    } else {
      /* Update the buffer's position. */
      map += h->masks[global].size;
    }
  }

  return map - reader->log.log.map;
}

/**
 * @brief Read the special flag of a particle (of any type) in the log file.
 *
 * @param reader The #logger_reader.
 * @param offset offset of the record to read.
 * @param mask (out) The mask of the record.
 * @param h_offset (out) Difference of offset with the next record.
 * @param data (out) The data of the flag.
 *
 * @return The special flag.
 */
__attribute__((always_inline)) INLINE enum logger_special_flags
logger_particle_read_special_flag(const struct logger_reader *reader,
                                  size_t offset, size_t *mask, size_t *h_offset,
                                  int *data) {

  /* Get a few pointers. */
  const struct header *h = &reader->log.header;
  void *map = reader->log.log.map;

  *mask = 0;
  *h_offset = 0;

  /* Read the record's mask. */
  map = logger_loader_io_read_mask(h, (char *)map + offset, mask, h_offset);

  /* Check that the mask is within the limits. */
  if (*mask > (unsigned int)(1 << h->masks_count)) {
    error_python("Found an unexpected mask %zi", *mask);
  }

  /* Check if it is not a time record. */
  if (*mask == h->timestamp_mask) {
    error_python("Cannot read a particle from timestep record.");
  }

  /* Read the special flag */
  uint32_t packed_data = 0;
  map = logger_loader_io_read_data(
      map, h->masks[logger_index_special_flags].size, &packed_data);

  return logger_unpack_flags_and_data(packed_data, data);
}
