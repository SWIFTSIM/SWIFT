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
 * @param all_fields The list of fields for the particle type.
 * @param all_fields_count The number of element in all_fields.
 * @param global_wanted The global index of the requested field.
 * @param mask (out) The mask of the record.
 * @param h_offset (out) Difference of offset with the next record.
 *
 * @return position after the record.
 */
__attribute__((always_inline)) INLINE size_t logger_particle_read_field(
    const struct logger_reader *reader, size_t offset, void *output,
    const struct field_information *all_fields, const int all_fields_count,
    const int global_wanted, size_t *mask, size_t *h_offset) {

  /* Get a few pointers. */
  const struct header *h = &reader->log.header;
  char *map = reader->log.log.map;

  *mask = 0;
  *h_offset = 0;

  /* Read the record's mask. */
  map = logger_loader_io_read_mask(h, map + offset, mask, h_offset);

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
  for (int i = 0; i < all_fields_count; i++) {
    const int global = all_fields[i].global_index;

    /* Is the mask present? */
    if (!(*mask & h->masks[global].mask)) {
      continue;
    }

    if (global_wanted == global) {
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
  char *map = reader->log.log.map;

  *mask = 0;
  *h_offset = 0;

  /* Read the record's mask. */
  map = logger_loader_io_read_mask(h, map + offset, mask, h_offset);

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

/**
 * @brief Interpolate a field of the particle at the given time.
 *
 * @param before Pointer to the #logger_field at a time < t.
 * @param after Pointer to the #logger_field at a time > t.
 * @param otuput Pointer to the output value.
 * @param time_before Time of field_before (< t).
 * @param time_after Time of field_after (> t).
 * @param time Requested time.
 * @param field The field to reconstruct.
 * @param params The simulation's #logger_parameters.
 */
void logger_particle_interpolate_field(
    const double time_before, const struct logger_field *restrict before,
    const double time_after, const struct logger_field *restrict after,
    void *restrict output, const double time,
    const struct field_information *field, enum part_type type,
    const struct logger_parameters *params) {

  /* Select the correct interpolation */
  switch (type) {

    /* Hydro */
    case swift_type_gas:
      switch (field->module) {
        case field_module_default:
          hydro_logger_interpolate_field(time_before, before, time_after, after,
                                         output, time, field->local_index,
                                         params);
          break;

        case field_module_chemistry:
          chemistry_logger_interpolate_field_part(
              time_before, before, time_after, after, output, time,
              field->local_index, params);
          break;

        default:
          error("Module not implemented");
      }
      break;

    /* Dark matter */
    case swift_type_dark_matter:
    case swift_type_dark_matter_background:
      if (field->module != field_module_default)
        error("Module not implemented");
      gravity_logger_interpolate_field(time_before, before, time_after, after,
                                       output, time, field->local_index,
                                       params);
      break;

    /* Stars */
    case swift_type_stars:
      switch (field->module) {
        case field_module_default:
          stars_logger_interpolate_field(time_before, before, time_after, after,
                                         output, time, field->local_index,
                                         params);
          break;

        case field_module_chemistry:
          chemistry_logger_interpolate_field_spart(
              time_before, before, time_after, after, output, time,
              field->local_index, params);
          break;

        case field_module_star_formation:
          star_formation_logger_interpolate_field(
              time_before, before, time_after, after, output, time,
              field->local_index, params);
          break;

        default:
          error("Module not implemented");
      }
      break;

      /* Default */
    default:
      error_python("Particle type not implemented");
  }
}
