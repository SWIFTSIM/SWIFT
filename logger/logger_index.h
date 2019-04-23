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

#include "logger_tools.h"

/**
 * @file logger_index.h
 * @brief This file deals with the index files.
 */
#ifndef __LOGGER_LOGGER_INDEX_H__
#define __LOGGER_LOGGER_INDEX_H__

struct logger_reader;


struct logger_index_data {
  /* The particle's id */
  long long id;

  /* The particle's offset */
  size_t offset;
};

/**
 * @brief This structure will contain the data related to
 *   the index file.
 *
 * It is initialized with #logger_index_init and freed with
 * #logger_index_free.
 *
 * The files are read with #logger_index_read_file and
 * can be freed with #logger_index_free_current_file.
 */
struct logger_index {
  /* The reader */
  struct logger_reader *reader;

  /* List of the time for each file */
  double *times;

  /* Number of files */
  int number_files;

  /* List of all the index filenames */
  char *basename;

  /* Particles' ids and offsets */
  struct logger_index_data *data;

  /* Number of particles */
  size_t total_number_particles;

  /* Number of particles per type */
  long long number_particles[swift_type_count];
};

void logger_index_init(struct logger_index *index, struct logger_reader *reader,
		       char *basename);
void logger_index_read_file(struct logger_index *index, int i);

void logger_index_free(struct logger_index *index);

void logger_index_free_current_file(struct logger_index *index);

#endif // __LOGGER_LOGGER_INDEX_H__
