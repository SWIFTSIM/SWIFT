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

/* Config parameters. */
#include "../config.h"

/* Some standard headers. */
#include <stdio.h>

/* This object's header. */
#include "logger_index.h"

#include "logger_loader_io.h"
#include "logger_tools.h"

/**
 * @brief Initialize the #logger_index by reading the index file.
 *
 * @param index The #logger_index.
 * @param reader The #logger_reader.
 * @param filename The filename.
 */
void logger_index_init(struct logger_index *index, struct logger_reader *reader,
                       char *filename) {

  /* Open file. */
  index->data = logger_loader_io_mmap_file(filename, &index->file_size,
                                           /* read_only */ 1);

  /* Read the double time. */
  size_t offset = 0;
  void *map = index->data + offset;
  map = logger_loader_io_read_data(map, sizeof(double), &index->time);

  /* Read the integer time. */
  map =
      logger_loader_io_read_data(map, sizeof(integertime_t), &index->int_time);

  /* Read the number of particles. */
  map = logger_loader_io_read_data(map, swift_type_count * sizeof(long long),
                                   &index->number_particles);

  /* Count total number of particles. */
  long long N = 0;
  for (int j = 0; j < swift_type_count; j++) {
    N += index->number_particles[j];
  }

  index->total_number_particles = N;
}

/**
 * @brief Free the memory.
 *
 * @param index The #logger_index.
 */
void logger_index_free(struct logger_index *index) {

  /* unmap file. */
  logger_loader_io_munmap_file(index->data, index->file_size);

  /* Set variables to default value. */
  index->total_number_particles = 0;

  for (int i = 0; i < swift_type_count; i++) {
    index->number_particles[i] = 0;
  }
}
