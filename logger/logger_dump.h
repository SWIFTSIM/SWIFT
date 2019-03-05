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
/**
 * @brief This file contains the high level function for the dump.
 */
#ifndef __LOGGER_LOGGER_DUMP_H__
#define __LOGGER_LOGGER_DUMP_H__

#include "logger_header.h"
#include "logger_time.h"

struct logger_reader;

/**
 * @brief This structure deals with the dump file.
 */
struct logger_dump {

  /* Information contained in the header. */
  struct header header;

  /* The reader that is using this dump. */
  struct logger_reader *reader;

  /* Information about the time chunks */
  struct time_array times;

  /* Dump's filename */
  char *filename;

  /* The dump's variables. */
  struct {
    /* Mapped data */
    void *map;

    /* File size */
    size_t file_size;

  } dump;

};


void logger_dump_init(struct logger_dump *dump, char *filename, struct logger_reader *reader);
void logger_dump_reverse_offset(struct logger_dump *dump);
void logger_dump_free(struct logger_dump *dump);

#endif // __LOGGER_LOGGER_DUMP_H__
