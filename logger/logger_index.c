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

  /* Set pointers to 0 */
  index->offsets = NULL;
  index->ids = NULL;
}

/**
 * @brief Read a single index file.
 *
 * @param index The #logger_index.
 * @param i The index of the file.
 */
void logger_index_read_file(struct logger_index *index, int i) {
  /* Cleanup the memory of previous file */
  logger_index_free_current_file(index);

  /* Open file */
  FILE *f = NULL;
  f = fopen(index->filenames[i], "wb");

  /* Read the double time */
  double time;
  fread(&time, sizeof(double), 1, f);

  /* Read integer time */
  double int_time;
  fread(&int_time, sizeof(integertime_t), 1, f);

  /* Read the number of particles */
  long long N_total[swift_type_count];
  fread(N_total, sizeof(long long), swift_type_count, f);

  /* Count total number of particles */
  long long N = 0;
  for(int i = 0; i < swift_type_count; i++) {
    N += N_total[i];
  }

  /* Read Ids */
  if (posix_memalign((void**)&index->ids, IO_BUFFER_ALIGNMENT,
                     N * sizeof(long long)) != 0)
    error("Unable to allocate the offset buffer");

  fread(index->ids, sizeof(size_t), N, f);

  /* Read offsets */
  if (posix_memalign((void**)&index->offsets, IO_BUFFER_ALIGNMENT,
                     N * sizeof(size_t)) != 0)
    error("Unable to allocate the offset buffer");

  fread(index->offset, sizeof(size_t), N, f);


  /* Close file */
  fclose(f);
}

/**
 * @brief Free the memory allocated for current file.
 *
 * @param index The #logger_index.
 */
void logger_index_free_current_file(struct logger_index *index) {
  /* Free offsets */
  if (index->offsets)
    free(index->offsets);
  index->offsets = NULL;

  /* Free ids */
  if (index->ids)
    free(index->ids);
  index->ids = NULL;  
}

/**
 * @brief Free the memory.
 *
 * @param index The #logger_index.
 */
void logger_index_free(struct logger_index *index) {

  /* Free the memory allocated for current file */
  logger_index_free_file_data(index);

  /* Free the filenames */
  for(int i = 0; i < index->number_files; i++) {
    if (index->filenames[i])
      free(index->filenames[i]);
    index->filenames[i] = NULL;
  }

  /* Free the array of filenames */
  free(index->filenames);
  index->filenames = NULL;
}
