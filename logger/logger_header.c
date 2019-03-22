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
#include "logger_header.h"

#include "logger_io.h"
#include "logger_tools.h"
#include "logger_dump.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Offset type names. */
const char *logger_offset_name[logger_offset_count] = {
  "Forward", "Backward", "Corrupted",
};



/**
 * @brief Print the properties of the header to stdout.
 *
 * @param h The #header.
 */
void header_print(const struct header *h) {
#ifdef SWIFT_DEBUG_CHECKS
  message("Debug checks enabled");
#endif
  message("Version:          %s", h->version);
  message("First Offset:     %lu", h->offset_first);
  message("Offset direction: %s", logger_offset_name[h->offset_direction]);
  message("Number masks:     %lu", h->number_mask);

  for (size_t i = 0; i < h->number_mask; i++) {
    message("\tMask:  %s", h->masks[i].name);
    message("\tValue: %u", h->masks[i].mask);
    message("\tSize:  %i", h->masks[i].size);
    message("");
  }
};

/**
 * @brief free allocated memory
 *
 * @param h The #header
 */
void header_free(struct header *h) { free(h->masks); };

/**
 * @brief Check if a field is present in the header.
 *
 * @param h The #header.
 * @param field name of the requested field.
 * @return Index of the field (-1 if not found).
 */
int header_get_field_index(const struct header *h, const char *field) {
  for (size_t i = 0; i < h->number_mask; i++) {
    if (strcmp(h->masks[i].name, field) == 0) {
      return i;
    }
  }

  return -1;
};

/**
 * @brief Inverse the offset direction
 *
 * @param h #header file structure
 * @param new_value The new value to write
 *
 */
void header_change_offset_direction(struct header *h, int new_value) {
  h->offset_direction = new_value;
  size_t offset = LOGGER_VERSION_SIZE;

  io_write_data(h->dump->dump.map, LOGGER_NUMBER_SIZE, &new_value, offset);
}

/**
 * @brief read the logger header
 *
 * @param h out: header
 * @param map file mapping
 */
void header_read(struct header *h, struct logger_dump *dump) {
  size_t offset = 0;
  void *map = dump->dump.map;

  /* Set pointer to dump */
  h->dump = dump;

  /* read version */
  offset = io_read_data(map, LOGGER_VERSION_SIZE, &h->version, offset);

  /* read offset direction */
  h->offset_direction = -1;
  offset = io_read_data(map, LOGGER_NUMBER_SIZE, &h->offset_direction, offset);

  if (!header_is_forward(h) && !header_is_backward(h) &&
      !header_is_corrupted(h))
    error("Wrong offset value in the header (%i)", h->offset_direction);

  /* read offset to first data */
  h->offset_first = 0;
  offset = io_read_data(map, LOGGER_OFFSET_SIZE, &h->offset_first, offset);

  /* read name size */
  h->name_length = 0;
  offset = io_read_data(map, LOGGER_NUMBER_SIZE, &h->name_length, offset);

  /* check if value defined in this file is large enough */
  if (STRING_SIZE < h->name_length) {
    error("Name too large in dump file");
  }

  /* read number of masks */
  h->number_mask = 0;
  offset = io_read_data(map, LOGGER_NUMBER_SIZE, &h->number_mask, offset);

  /* allocate memory */
  h->masks = malloc(sizeof(struct mask_data) * h->number_mask);

  /* loop over all masks */
  for (size_t i = 0; i < h->number_mask; i++) {
    /* read mask name */
    offset = io_read_data(map, h->name_length, h->masks[i].name, offset);

    /* get mask value */
    h->masks[i].mask = 1 << i;

    /* read mask data size */
    h->masks[i].size = 0;
    offset = io_read_data(map, LOGGER_NUMBER_SIZE, &h->masks[i].size, offset);
  }

  if (offset != h->offset_first) {
    header_print(h);
    error("Wrong header size (in header %li, current %li)", h->offset_first,
          offset);
  }
};

/**
 * @brief count number of bits in a given mask
 *
 * @param h #header file structure
 * @param mask mask to compute
 *
 * @return number of bits in mask
 */
size_t header_get_mask_size(const struct header *h, const size_t mask) {
  size_t count = 0;
  for (size_t i = 0; i < h->number_mask; i++) {
    if (mask & h->masks[i].mask) {
      count += h->masks[i].size;
    }
  }
  return count;
}
