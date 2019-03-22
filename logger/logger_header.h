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
#ifndef __LOGGER_LOGGER_HEADER_H__
#define __LOGGER_LOGGER_HEADER_H__

#include "logger_tools.h"

#include <stdio.h>
#include <stdlib.h>

#define LOGGER_VERSION_SIZE 20
#define LOGGER_NUMBER_SIZE 4
#define LOGGER_OFFSET_SIZE 7
#define LOGGER_MASK_SIZE 1

enum logger_offset_direction {
  logger_offset_backward = 0,
  logger_offset_forward,
  logger_offset_corrupted,
  /* Number of offset type */
  logger_offset_count,
};

/**
 * @brief Names of the offset directions.
 */
extern const char *logger_offset_name[];


struct logger_dump;

/**
 * @brief This structure contains everything from the file header.
 *
 * This structure is initialized by @header_read and need to be freed
 * with @header_free.
 *
 * The information contained by the header can be easily access with
 * the functions @header_get_mask_size and @header_get_field_index.
 *
 * The only function that modify the file is @header_change_offset_direction.
 */
struct header {
  /* Logger version. */
  char version[STRING_SIZE];

  /* Offset of the first header. */
  size_t offset_first;

  /* Number of bytes for names. */
  size_t name_length;

  /* Number of masks. */
  size_t number_mask;

  /* List of masks. */
  struct mask_data *masks;

  /* Direction of the offset in the chunks. */
  int offset_direction;

  /* The corresponding dump */
  struct logger_dump *dump;
};

void header_print(const struct header *h);
void header_free(struct header *h);
int header_get_field_index(const struct header *h, const char *field);
void header_read(struct header *h, struct logger_dump *dump);
size_t header_get_mask_size(const struct header *h, const size_t mask);
void header_change_offset_direction(struct header *h, int new_value);


/**
 * @brief Check if the offset are forward.
 * @param h The #header.
 */
__attribute__((always_inline)) INLINE static int header_is_forward(
    const struct header *h) {
  return h->offset_direction == logger_offset_forward;
}

/**
 * @brief Check if the offset are backward.
 * @param h The #header.
 */
__attribute__((always_inline)) INLINE static int header_is_backward(
    const struct header *h) {
  return h->offset_direction == logger_offset_backward;
}

/**
 * @brief Check if the offset are corrupted.
 * @param h The #header.
 */
__attribute__((always_inline)) INLINE static int header_is_corrupted(
    const struct header *h) {
  return h->offset_direction == logger_offset_corrupted;
}

#endif  // __LOGGER_LOGGER_HEADER_H__
