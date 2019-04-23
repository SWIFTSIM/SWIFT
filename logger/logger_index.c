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

#include "logger_tools.h"

/**
 * @brief Initialize the #logger_index.
 *
 * @param index The #logger_index.
 * @param reader The #logger_reader.
 * @param basename The basename of the index files.
 */
void logger_index_init(struct logger_index *index, struct logger_reader *reader,
		       char *basename) {

  /* Set pointers to 0. */
  index->data = NULL;

  /* Set variables default value. */
  index->total_number_particles = 0;

  for(int i = 0; i < swift_type_count; i++) {
    index->number_particles[i] = 0;
  }
}

/**
 * @brief Read a single index file.
 *
 * @param index The #logger_index.
 * @param i The index of the file.
 */
void logger_index_read_file(struct logger_index *index, int i) {
  /* Cleanup the memory of previous file. */
  logger_index_free_current_file(index);

  /* Open file. */
  FILE *f = NULL;
  char name[200];
  sprintf(name, "%s_%04i.index", index->basename, i);
  f = fopen(name, "rb");

  /* Read the double time. */
  double time;
  fread(&time, sizeof(double), 1, f);

  /* Read the integer time. */
  double int_time;
  fread(&int_time, sizeof(integertime_t), 1, f);

  /* Read the number of particles. */
  fread(index->number_particles, sizeof(long long), swift_type_count, f);

  /* Count total number of particles. */
  long long N = 0;
  for(int j = 0; j < swift_type_count; j++) {
    N += index->number_particles[j];
  }

  index->total_number_particles = N;


  /* Read the particles ids. */
  if (posix_memalign((void**)&index->data, IO_BUFFER_ALIGNMENT,
                     N * sizeof(struct logger_index_data)) != 0)
    error("Unable to allocate index data buffer");

  fread(index->data, sizeof(struct logger_index_data), N, f);

  /* Close the file. */
  fclose(f);
}

/**
 * @brief Free the memory allocated for current file.
 *
 * @param index The #logger_index.
 */
void logger_index_free_current_file(struct logger_index *index) {
  /* Free the index data */
  if (index->data)
    free(index->data);
  index->data = NULL;

  /* Set variables to default value. */
  index->total_number_particles = 0;

  for(int i = 0; i < swift_type_count; i++) {
    index->number_particles[i] = 0;
  }
}

/**
 * @brief Free the memory.
 *
 * @param index The #logger_index.
 */
void logger_index_free(struct logger_index *index) {

  /* Free the memory allocated for current file. */
  logger_index_free_current_file(index);
}
