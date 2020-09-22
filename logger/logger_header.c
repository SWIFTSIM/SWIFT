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

#include "logger_loader_io.h"
#include "logger_logfile.h"
#include "logger_reader.h"
#include "logger_tools.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* Name of each offset direction. */
const char *logger_offset_name[logger_offset_count] = {
    "Forward",
    "Backward",
    "Corrupted",
};

/**
 * @brief Print the properties of the header to stdout.
 *
 * @param h The #header.
 */
void header_print(const struct header *h) {
#ifdef SWIFT_DEBUG_CHECKS
  message("Debug checks enabled.");
#endif
  message("First Offset:     %lu.", h->offset_first_record);
  message("Offset direction: %s.", logger_offset_name[h->offset_direction]);
  message("Number masks:     %i.", h->masks_count);

  for (int i = 0; i < h->masks_count; i++) {
    message("  Mask:  %s.", h->masks[i].name);
    message("  Value: %u.", h->masks[i].mask);
    message("  Size:  %i.", h->masks[i].size);
    message("");
  }
};

/**
 * @brief free the allocated memory.
 *
 * @param h The #header.
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
  for (int i = 0; i < h->masks_count; i++) {
    if (strcmp(h->masks[i].name, field) == 0) {
      return i;
    }
  }

  return -1;
};

/**
 * @brief Update the offset direction in the structure and
 * write it to the logfile.
 *
 * @param h #header file structure.
 * @param new_value The new value to write.
 *
 */
void header_change_offset_direction(struct header *h,
                                    enum logger_offset_direction new_value) {
  h->offset_direction = new_value;
  /* Skip file format and version numbers. */
  size_t offset = LOGGER_VERSION_SIZE + 2 * sizeof(int);

  logger_loader_io_write_data((char *)h->log->log.map + offset,
                              sizeof(unsigned int), &new_value);
}

/**
 * @brief read the logger header.
 *
 * @param h out: The #header.
 * @param log The #logger_logfile.
 */
void header_read(struct header *h, struct logger_logfile *log) {
  void *map = log->log.map;

  /* Set pointer to log. */
  h->log = log;

  /* read the file format. */
  char file_format[STRING_SIZE];
  map = logger_loader_io_read_data(map, LOGGER_VERSION_SIZE, &file_format);
  if (strcmp(file_format, "SWIFT_LOGGER"))
    error_python("Wrong file format (%s).", file_format);

  /* Read the major version number. */
  map = logger_loader_io_read_data(map, sizeof(int), &h->major_version);

  /* Read the minor version number. */
  map = logger_loader_io_read_data(map, sizeof(int), &h->minor_version);

  struct logger_reader *reader = log->reader;
  if (&reader->log != log) error_python("Wrong link to the reader.");

  if (reader->verbose > 0)
    message("File version %i.%i.", h->major_version, h->minor_version);

  /* Read the offset directions. */
  map = logger_loader_io_read_data(map, sizeof(int), &h->offset_direction);

  if (!header_is_forward(h) && !header_is_backward(h) &&
      !header_is_corrupted(h))
    error_python("Wrong offset value in the header (%i).", h->offset_direction);

  /* Read offset to first record. */
  h->offset_first_record = 0;
  map = logger_loader_io_read_data(map, LOGGER_OFFSET_SIZE,
                                   &h->offset_first_record);

  /* Read the size of the strings. */
  h->string_length = 0;
  map =
      logger_loader_io_read_data(map, sizeof(unsigned int), &h->string_length);

  /* Check if value defined in this file is large enough. */
  if (STRING_SIZE < h->string_length) {
    error_python("Name too large in log file %i.", h->string_length);
  }

  /* Read the number of masks. */
  map = logger_loader_io_read_data(map, sizeof(unsigned int), &h->masks_count);

  /* Allocate the masks memory. */
  h->masks = malloc(sizeof(struct mask_data) * h->masks_count);
  if (h->masks == NULL) {
    error_python("Failed to allocate the memory for the masks.");
  }

  /* Loop over all masks. */
  h->timestamp_mask = 0;
  for (int i = 0; i < h->masks_count; i++) {
    /* Read the mask name. */
    map = logger_loader_io_read_data(map, h->string_length, h->masks[i].name);

    /* Set the mask value. */
    h->masks[i].mask = 1 << i;

    /* Read the mask data size. */
    map = logger_loader_io_read_data(map, sizeof(unsigned int),
                                     &h->masks[i].size);

    /* Keep the timestamp mask in memory */
    if (strcmp(h->masks[i].name, "Timestamp") == 0) {
      h->timestamp_mask = h->masks[i].mask;
    }
  }

  /* Check that the timestamp mask exists */
  if (h->timestamp_mask == 0) {
    error_python("Unable to find the timestamp mask.");
  }

  /* Check the logfile header's size. */
  if (map != (char *)log->log.map + h->offset_first_record) {
    size_t offset = (char *)map - (char *)log->log.map;
    error_python("Wrong header size (in header %zi, current %zi).",
                 h->offset_first_record, offset);
  }

  /* Set the first and second derivatives as non existent. */
  for (int i = 0; i < h->masks_count; i++) {
    h->masks[i].reader.first_deriv = -1;
    h->masks[i].reader.second_deriv = -1;
  }

  /* Hydro */
  /* Set the link between local and global */
  for (int j = 0; j < hydro_logger_field_count; j++) {
    hydro_logger_local_to_global[j] = -1;
    for (int i = 0; i < h->masks_count; i++) {
      if (strcmp(h->masks[i].name, hydro_logger_field_names[j]) == 0) {
        hydro_logger_local_to_global[j] = i;
        break;
      }
    }

    /* Check if everything is fine. */
    const int index = hydro_logger_local_to_global[j];
    if (index == -1) {
      error("Field %s in hydro is not set", hydro_logger_field_names[j]);
    }
    if (h->masks[index].size != hydro_logger_field_size[j]) {
      error("Field %s in hydro does not have the correct size",
            hydro_logger_field_names[j]);
    }
  }

  hydro_logger_reader_link_derivatives(h);

  /* Gravity */
  /* Set the link between local and global */
  for (int j = 0; j < gravity_logger_field_count; j++) {
    gravity_logger_local_to_global[j] = -1;
    for (int i = 0; i < h->masks_count; i++) {
      if (strcmp(h->masks[i].name, gravity_logger_field_names[j]) == 0) {
        gravity_logger_local_to_global[j] = i;
        break;
      }
    }

    /* Check if everything is fine. */
    const int index = gravity_logger_local_to_global[j];
    if (index == -1) {
      error("Field %s in gravity is not set", gravity_logger_field_names[j]);
    }
    if (h->masks[index].size != gravity_logger_field_size[j]) {
      error("Field %s in gravity does not have the correct size",
            gravity_logger_field_names[j]);
    }
  }

  gravity_logger_reader_link_derivatives(h);

  /* Stars */
  /* Set the link between local and global */
  for (int j = 0; j < stars_logger_field_count; j++) {
    stars_logger_local_to_global[j] = -1;
    for (int i = 0; i < h->masks_count; i++) {
      if (strcmp(h->masks[i].name, stars_logger_field_names[j]) == 0) {
        stars_logger_local_to_global[j] = i;
        break;
      }
    }

    /* Check if everything is fine. */
    const int index = stars_logger_local_to_global[j];
    if (index == -1) {
      error("Field %s in stars is not set", stars_logger_field_names[j]);
    }
    if (h->masks[index].size != stars_logger_field_size[j]) {
      error("Field %s in stars does not have the correct size.",
            stars_logger_field_names[j]);
    }
  }
  stars_logger_reader_link_derivatives(h);
};

/**
 * @brief Count number of bits in a given mask (without the record header).
 *
 * @param h #header file structure.
 * @param mask Mask to compute.
 *
 * @return number of bits in mask.
 */
size_t header_get_record_size_from_mask(const struct header *h,
                                        const size_t mask) {
  size_t count = 0;
  /* Loop over each masks. */
  for (int i = 0; i < h->masks_count; i++) {
    if (mask & h->masks[i].mask) {
      count += h->masks[i].size;
    }
  }
  return count;
}
