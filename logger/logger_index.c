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

/* Include the corresponding header */
#include "logger_index.h"

/* Include the standard headers */
#include <errno.h>
#include <fcntl.h>

/* Include local headers */
#include "logger_loader_io.h"
#include "logger_reader.h"
#include "quick_sort.h"

/* Define the place and size of each element in the header. */
/* The time in double precision. */
#define logger_index_time_offset 0
#define logger_index_time_size sizeof(double)
/* The time on the integer timeline. */
#define logger_index_integer_time_offset \
  logger_index_time_offset + logger_index_time_size
#define logger_index_integer_time_size sizeof(integertime_t)
/* The number of particle (for each type). */
#define logger_index_npart_offset \
  logger_index_integer_time_offset + logger_index_integer_time_size
#define logger_index_npart_size sizeof(uint64_t) * swift_type_count
/* The flag used for checking if the file is sorted or not. */
#define logger_index_is_sorted_offset \
  logger_index_npart_offset + logger_index_npart_size
#define logger_index_is_sorted_size sizeof(char)
/* The array containing the offset and ids. Rounding to a multiplier of 8 in
 * order to align the memory */
#define logger_index_data_offset \
  ((logger_index_is_sorted_offset + logger_index_is_sorted_size + 7) & ~7)

/**
 * @brief Read the index file header.
 *
 * @param index The #logger_index.
 * @param filename The filename of the index file.
 */
void logger_index_read_header(struct logger_index *index,
                              const char *filename) {

  /* Open the file */
  if (index->reader->verbose > 1) {
    message("Reading %s", filename);
  }
  logger_index_map_file(index, filename, 0);

  /* Read times */
  memcpy(&index->time, (char *)index->index.map + logger_index_time_offset,
         logger_index_time_size);
  memcpy(&index->integer_time,
         (char *)index->index.map + logger_index_integer_time_offset,
         logger_index_integer_time_size);

  /* Read the number of particles. */
  memcpy(index->nparts, (char *)index->index.map + logger_index_npart_offset,
         logger_index_npart_size);

  /* Read whether the file is sorted. */
  memcpy(&index->is_sorted,
         (char *)index->index.map + logger_index_is_sorted_offset,
         logger_index_is_sorted_size);

  /* Compute the position of the history of new particles. */
  size_t count = 0;
  for (int i = 0; i < swift_type_count; i++) {
    count += index->nparts[i];
  }
  count *= sizeof(struct index_data);

  /* Read the number of new particles. */
  memcpy(&index->nparts_created,
         (char *)index->index.map + logger_index_data_offset + count,
         logger_index_npart_size);

  /* Compute the position of the history of particles removed. */
  size_t count_new = 0;
  for (int i = 0; i < swift_type_count; i++) {
    count_new += index->nparts_created[i];
  }
  count_new *= sizeof(struct index_data);

  /* Read the number of particles removed. */
  memcpy(&index->nparts_removed,
         (char *)index->index.map + logger_index_data_offset + count +
             logger_index_npart_size + count_new,
         logger_index_npart_size);

  /* Close the file */
  logger_index_free(index);
}

/**
 * @brief Write that the file is sorted.
 *
 * WARNING The file must be mapped.
 *
 * @param index The #logger_index.
 */
void logger_index_write_sorted(struct logger_index *index) {
  /* Set the value */
  const char is_sorted = 1;

  /* Write the value */
  memcpy((char *)index->index.map + logger_index_is_sorted_offset, &is_sorted,
         logger_index_is_sorted_size);
}

/**
 * @brief Get the #index_data of a particle type.
 *
 * The #index_data contains the ids and offset of the particle.
 * This function should be used when looking for the last known particles.
 *
 * @param index The #logger_index.
 * @param type The particle type.
 */
struct index_data *logger_index_get_data(struct logger_index *index, int type) {
  /* Get the offset of particles of this type by skipping entries of lower
   * types. */
  size_t count = 0;
  for (int i = 0; i < type; i++) {
    count += index->nparts[i];
  }
  const size_t type_offset = count * sizeof(struct index_data);

  const char *ret =
      (char *)index->index.map + logger_index_data_offset + type_offset;
  return (struct index_data *)ret;
}

/**
 * @brief Get the #index_data for the particles created.
 *
 * The #index_data contains the ids and offset of the particle.
 * This function should be used when looking for the particles created.
 *
 * @param index The #logger_index.
 * @param type The particle type.
 */
struct index_data *logger_index_get_created_history(struct logger_index *index,
                                                    int type) {

  /* Get the position after the previous array. */
  char *ret = (char *)logger_index_get_data(index, swift_type_count);

  /* Get the offset of particles of this type by skipping entries of lower
   * types. */
  size_t count = 0;
  for (int i = 0; i < type; i++) {
    count += index->nparts_created[i];
  }
  const size_t type_offset = count * sizeof(struct index_data);

  ret += type_offset + logger_index_npart_size;
  return (struct index_data *)ret;
}

/**
 * @brief Get the #index_data for the particles removed.
 *
 * The #index_data contains the ids and offset of the particle.
 * This function should be used when looking for the particles removed.
 *
 * @param index The #logger_index.
 * @param type The particle type.
 */
struct index_data *logger_index_get_removed_history(struct logger_index *index,
                                                    int type) {

  /* Get the position after the previous array. */
  char *ret = (char *)logger_index_get_created_history(index, swift_type_count);

  /* Get the offset of particles of this type by skipping entries of lower
   * types. */
  size_t count = 0;
  for (int i = 0; i < type; i++) {
    count += index->nparts_removed[i];
  }
  const size_t type_offset = count * sizeof(struct index_data);

  ret += type_offset + logger_index_npart_size;
  return (struct index_data *)ret;
}

/**
 * @brief Check if a time array is written at the end.
 *
 * @param index The #logger_index.
 */
int logger_index_contains_time_array(struct logger_index *index) {
  /* Only the first index file should have a time array */
  if (index->time != 0.) {
    error_python("Only the first index file can have a time array.");
  }

  /* Get the file size */
  const size_t file_size = index->index.mmap_size;

  /* Get the position after the previous array. */
  char *ret = (char *)logger_index_get_removed_history(index, swift_type_count);
  const size_t before_time_array = ret - (char *)index->index.map;

  return file_size != before_time_array;
}

/**
 * @brief Map the file and if required sort it.
 *
 * @param index The #logger_index.
 * @param filename The index filename.
 * @param sorted Does the file needs to be sorted?
 */
void logger_index_map_file(struct logger_index *index, const char *filename,
                           int sorted) {
  /* Un-map previous file */
  if (index->index.map != NULL) {
    error_python("Trying to remap.");
  }

  /* Check if need to sort the file */
  if (sorted && !index->is_sorted) {
    if (index->reader->verbose > 1) {
      message("Sorting the index file.");
    }
    /* Map the index file */
    logger_loader_io_mmap_file(&index->index, filename,
                               /* read_only */ 0);
    /* Sort the file */
    for (int i = 0; i < swift_type_count; i++) {
      if (index->nparts[i] != 0) {
        struct index_data *data = logger_index_get_data(index, i);
        quick_sort(data, index->nparts[i]);
      }
    }

    /* Write that the file is sorted */
    logger_index_write_sorted(index);

    /* Free the index file before opening it again in read only */
    logger_index_free(index);

    if (index->reader->verbose > 1) {
      message("Sorting done.");
    }
  }

  /* Map the index file */
  logger_loader_io_mmap_file(&index->index, filename,
                             /* read_only */ 1);
}

/**
 * @brief Cleanup the memory of a logger_index
 *
 * @param index The #logger_index.
 */
void logger_index_free(struct logger_index *index) {
  if (index->index.map == NULL) {
    error_python("Trying to unmap an unexisting map");
  }
  logger_loader_io_munmap_file(&index->index);
}

/**
 * @brief Get the offset of a given particle
 *
 * @param index The #logger_index.
 * @param id The ID of the particle.
 * @param type The type of the particle.
 *
 * @return The offset of the particle or 0 if not found.
 */
size_t logger_index_get_particle_offset(struct logger_index *index,
                                        long long id, int type) {
  /* Define a few variables */
  const struct index_data *data = logger_index_get_data(index, type);
  size_t left = 0;
  size_t right = index->nparts[type] - 1;

  /* Search for the value (binary search) */
  while (left <= right) {
    size_t m = (left + right) / 2;
    if (data[m].id < id) {
      left = m + 1;
    } else if (data[m].id > id) {
      right = m - 1;
    } else {
      return data[m].offset;
    }
  }

  return 0;
}

/**
 * @brief Initialize the #logger_index.
 *
 * @param index The #logger_index.
 * @param reader The #logger_reader.
 */
void logger_index_init(struct logger_index *index,
                       struct logger_reader *reader) {
  /* Set the mapped file to NULL */
  index->index.map = NULL;

  /* Set the pointer to the reader */
  index->reader = reader;

  /* Set the time to its default value */
  index->time = -1;
}
