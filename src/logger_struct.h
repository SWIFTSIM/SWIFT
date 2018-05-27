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

#define LOGGER_STRING_LENGTH 200
#include "dump.h"

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
