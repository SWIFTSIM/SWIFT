/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
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
#ifndef SWIFT_LOGGER_STRUCT_H
#define SWIFT_LOGGER_STRUCT_H

#ifdef WITH_LOGGER

#include <limits.h>

/* Local Includes */
#include "dump.h"

#define LOGGER_STRING_LENGTH 200

/* parameters of the logger */
struct logger_parameters {
  /* size of a label in bytes */
  size_t label_size;

  /* size of an offset in bytes */
  size_t offset_size;

  /* size of a mask in bytes */
  size_t mask_size;

  /* size of a number in bytes */
  size_t number_size;

  /* size of a data type in bytes */
  size_t data_type_size;
  
  /* number of different mask */
  size_t nber_mask;

  /* value of each masks */
  size_t *masks;

  /* data size of each mask */
  size_t *masks_data_size;
  
  /* label of each mask */
  char *masks_name;

};


/* structure containing global data */
struct logger {
  /* Number of particle steps between dumping a chunk of data */
  short int delta_step;

  /* Logger basename */
  char base_name[LOGGER_STRING_LENGTH];  

  /* File name of the dump file */
  struct dump *dump;

  /* timestamp offset for logger*/
  size_t timestamp_offset;

  /* size of the buffer */
  size_t buffer_size;

  /* scaling factor when buffer is too small */
  float buffer_scale;

  /* logger parameters */
  struct logger_parameters *params;

} SWIFT_STRUCT_ALIGN;

/* required structure for each particle type */
struct logger_part_data {
  /* Number of particle updates since last output */
  short int last_output;

  /* offset of last particle log entry */
  size_t last_offset;
};

__attribute__((always_inline)) INLINE static  void logger_part_data_init(
    struct logger_part_data *logger ) {
  logger->last_offset = 0;
  logger->last_output = SHRT_MAX;
}

#endif // WITH_LOGGER

#endif // SWIFT_LOGGER_STRUCT_H
