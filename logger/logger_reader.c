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

/* Include corresponding header */
#include "logger_reader.h"

/* Include standard library */
#include <sys/sysinfo.h>
#include <unistd.h>

/* Include local headers */
#include "threadpool.h"

#define nr_threads get_nprocs()

/**
 * @brief Initialize the reader.
 *
 * @param reader The #logger_reader.
 * @param basename The basename of the logger files.
 * @param verbose The verbose level.
 */
void logger_reader_init(struct logger_reader *reader, const char *basename,
                        int verbose) {
  if (verbose > 1) message("Initializing the reader.");

  /* Set the variable to the default values */
  reader->time.time = -1.;
  reader->time.int_time = 0;
  reader->time.time_offset = 0;

  /* Copy the base name */
  strcpy(reader->basename, basename);

  /* Initialize the reader variables. */
  reader->verbose = verbose;

  /* Generate the logfile filename */
  char logfile_name[STRING_SIZE];
  sprintf(logfile_name, "%s.dump", basename);

  /* Initialize the log file. */
  logger_logfile_init_from_file(&reader->log, logfile_name, reader,
                                /* only_header */ 0);

  /* Initialize the index files */
  logger_reader_init_index(reader);

  if (verbose > 1) message("Initialization done.");
}

/**
 * @brief Initialize the index part of the reader.
 *
 * @param reader The #logger_reader.
 */
void logger_reader_init_index(struct logger_reader *reader) {
  /* Initialize the logger_index */
  logger_index_init(&reader->index.index, reader);

  /* Count the number of files */
  int count = 0;
  while (1) {
    char filename[STRING_SIZE + 50];
    sprintf(filename, "%s_%04i.index", reader->basename, count);

    /* Check if file exists */
    if (access(filename, F_OK) != -1) {
      count++;
    } else {
      break;
    }
  }

  reader->index.n_files = count;

  /* Initialize the arrays */
  if ((reader->index.times = (double *)malloc(count * sizeof(double))) ==
      NULL) {
    error("Failed to allocate the list of times");
  }
  if ((reader->index.int_times =
           (integertime_t *)malloc(count * sizeof(integertime_t))) == NULL) {
    error("Failed to allocate the list of times");
  }

  /* Get the information contained in the headers */
  for (int i = 0; i < reader->index.n_files; i++) {
    char filename[STRING_SIZE + 50];
    sprintf(filename, "%s_%04i.index", reader->basename, i);

    /* Read the header */
    logger_index_read_header(&reader->index.index, filename);

    /* Save the required information */
    reader->index.times[i] = reader->index.index.time;
    reader->index.int_times[i] = reader->index.index.integer_time;
  }
}

/**
 * @brief Free the reader.
 *
 * @param reader The #logger_reader.
 */
void logger_reader_free(struct logger_reader *reader) {
  /* Free the log. */
  logger_logfile_free(&reader->log);

  if (reader->time.time != -1.) {
    logger_index_free(&reader->index.index);
  }

  free(reader->index.int_times);
  free(reader->index.times);
}

/**
 * @brief Set the reader to a given time and read the correct index file.
 *
 * @param reader The #logger_reader.
 * @param time The requested time.
 */
void logger_reader_set_time(struct logger_reader *reader, double time) {
  /* Set the time */
  reader->time.time = time;

  /* Find the correct index */
  unsigned int left = 0;
  unsigned int right = reader->index.n_files - 1;

  while (left != right) {
    /* Do a ceil - division */
    unsigned int m = (left + right + 1) / 2;
    if (reader->index.times[m] > time) {
      right = m - 1;
    } else {
      left = m;
    }
  }

  /* Generate the filename */
  char filename[STRING_SIZE + 50];
  sprintf(filename, "%s_%04u.index", reader->basename, left);

  /* Check if the file is already mapped */
  if (reader->index.index.index.map != NULL) {
    logger_index_free(&reader->index.index);
  }

  /* Read the file */
  logger_index_read_header(&reader->index.index, filename);
  logger_index_map_file(&reader->index.index, filename, /* sorted */ 1);

  /* Get the offset of the time chunk */
  size_t ind = time_array_get_index_from_time(&reader->log.times, time);
  /* ind == 0 and ind == 1 are the same time, but when reading we need
     data before the initial time. */
  if (ind == 0) {
    ind = 1;
  }

  /* Check if we requested exactly a time step  */
  if (reader->log.times.records[ind].time != time) {
    /* In order to interpolate, we need to be above and not below the time */
    ind += 1;
  }

  /* Save the values */
  reader->time.index = ind;
  reader->time.int_time = reader->log.times.records[ind].int_time;
  reader->time.time_offset = reader->log.times.records[ind].offset;
}

/**
 * @brief Provides the number of particle (per type) from the index file.
 *
 * @param reader The #logger_reader.
 * @param n_type (output) The number of particle type possible.
 *
 * @return For each type possible, the number of particle.
 */
const uint64_t *logger_reader_get_number_particles(struct logger_reader *reader,
                                                   int *n_type) {
  *n_type = swift_type_count;
  return reader->index.index.nparts;
}

/**
 * @brief Read a single field for a single part from the logger and interpolate
 it if needed.
 * This function will jump from the last known full record until the requested
 time and then
 * read the last record containing the requested field.
 * If an interpolation is required, the function will continue jumping until
 finding a record
 * with the requested field. If the first or second derivative is present in the
 last record before /
 * first record after the requested time, it will use it for the interpolation.
 *
 * @param reader The #logger_reader.
 * @param time The requested time.
 * @param offset_time The offset of the corresponding time record.
 * @param interp_type The type of interpolation requested.
 * @param offset_last_full_record The offset of this particle last record
 containing all the fields.
 * @param field The fields wanted in local index.
 * @param first_deriv_wanted The field that corresponds to the first derivative
 of fields (-1 if none)
 * @param second_deriv_wanted The field that corresponds to the second
 derivative of field (-1 if none)
 * @param output Pointers to the output array
 */
void logger_reader_read_field(struct logger_reader *reader, double time,
                              size_t offset_time,
                              enum logger_reader_type interp_type,
                              const size_t offset_last_full_record,
                              const int field, const int first,
                              const int second, void *output,
                              enum part_type type) {

  const struct header *h = &reader->log.header;
  size_t offset = offset_last_full_record;

  /* Get the particle type variables. */
  int *local_to_global = NULL;
  int field_count = 0;
  switch (type) {
    case swift_type_gas:
      local_to_global = hydro_logger_local_to_global;
      field_count = hydro_logger_field_count;
      break;
    case swift_type_dark_matter:
      local_to_global = gravity_logger_local_to_global;
      field_count = gravity_logger_field_count;
      break;
    case swift_type_stars:
      local_to_global = stars_logger_local_to_global;
      field_count = stars_logger_field_count;
      break;
    default:
      error("Particle type not implemented");
  }

  /* Get the masks. */
  struct mask_data *mask_field = &h->masks[local_to_global[field]];
  struct mask_data *mask_first = NULL;
  if (first >= 0) {
    mask_first = &h->masks[local_to_global[first]];
  }
  struct mask_data *mask_second = NULL;
  if (second >= 0) {
    mask_second = &h->masks[local_to_global[second]];
  }

  /* Offset of the field read before the requested time. */
  size_t offset_before = offset_last_full_record;

  /* Record's header information */
  size_t mask, h_offset;

  /* Find the data for the previous record.
     As we start from a full record,
     no need to check if the field is found.
  */
  while (offset < offset_time) {
    /* Read the particle. */
    logger_loader_io_read_mask(h, reader->log.log.map + offset, &mask,
                               &h_offset);

    /* Check if the field is present */
    if (mask & mask_field->mask) {
      offset_before = offset;
    }

    /* Go to the next record. */
    offset += h_offset;
  }

  /* Read the field */
  logger_particle_read_field(reader, offset_before, output, local_to_global,
                             field_count, field, &mask, &h_offset);

  /* Deal with the first derivative. */
  const size_t size_first = mask_first == NULL ? 0 : mask_first->size;
  char first_deriv[size_first];

  int first_found = mask_first != NULL && mask & mask_first->mask;
  if (first_found) {
    /* Read the first derivative */
    logger_particle_read_field(reader, offset_before, first_deriv,
                               local_to_global, field_count, first, &mask,
                               &h_offset);
  }

  /* Deal with the second derivative. */
  const size_t size_second = mask_second == NULL ? 0 : mask_second->size;
  char second_deriv[size_second];

  int second_found = mask_second != NULL && mask & mask_second->mask;
  if (second_found) {
    /* Read the first derivative */
    logger_particle_read_field(reader, offset_before, second_deriv,
                               local_to_global, field_count, second, &mask,
                               &h_offset);
  }

  /* Get the time. */
  // TODO reduce search interval
  double time_before = time_array_get_time(&reader->log.times, offset_before);

  /* When interpolating, we need to get the next record after
     the requested time.
  */
  if (interp_type != logger_reader_const) {

    /* Loop over the records until having all the fields. */
    while (1) {

      /* Read the particle. */
      logger_loader_io_read_mask(h, reader->log.log.map + offset, &mask,
                                 &h_offset);

      /* Do we have the field? */
      if (mask & mask_field->mask) {
        break;
      }

      /* Check if we can still move forward */
      if (h_offset == 0) {
        error("There is no record after the current one");
      }

      /* Go to the next record. */
      offset += h_offset;
    }

    /* Output after the requested time. */
    char output_after[mask_field->size];

    /* Read the field */
    logger_particle_read_field(reader, offset, output_after, local_to_global,
                               field_count, field, &mask, &h_offset);

    /* Deal with the first derivative. */
    char first_deriv_after[size_first];

    /* Did we find the derivative before and in this record? */
    first_found = mask_first != NULL && first_found && mask & mask_first->mask;
    if (first_found) {
      /* Read the first derivative */
      logger_particle_read_field(reader, offset, first_deriv_after,
                                 local_to_global, field_count, first, &mask,
                                 &h_offset);
    }

    /* Deal with the second derivative. */
    char second_deriv_after[size_second];

    /* Did we find the derivative before and in this record? */
    second_found =
        mask_second != NULL && second_found && mask & mask_second->mask;
    if (second_found) {
      /* Read the second derivative */
      logger_particle_read_field(reader, offset, second_deriv_after,
                                 local_to_global, field_count, second, &mask,
                                 &h_offset);
    }

    /* Get the time. */
    // TODO reduce search interval
    double time_after = time_array_get_time(&reader->log.times, offset);

    /* Deal with the derivatives */
    struct logger_field before;
    struct logger_field after;
    before.field = output;
    before.first_deriv = first_found ? first_deriv : NULL;
    before.second_deriv = second_found ? second_deriv : NULL;
    after.field = output_after;
    after.first_deriv = first_found ? first_deriv_after : NULL;
    after.second_deriv = second_found ? second_deriv_after : NULL;

    /* Interpolate the data. */
    switch (type) {
      case swift_type_gas:
        hydro_logger_interpolate_field(time_before, &before, time_after, &after,
                                       output, time, field);
        break;
      case swift_type_dark_matter:
        gravity_logger_interpolate_field(time_before, &before, time_after,
                                         &after, output, time, field);
        break;
      case swift_type_stars:
        stars_logger_interpolate_field(time_before, &before, time_after, &after,
                                       output, time, field);
        break;
      default:
        error("Particle type not implemented");
    }
  }
}

/**
 * @brief Convert the fields from global indexes to local.
 *
 * @param reader The #logger_reader.
 * @param global_fields_wanted The fields to sort.
 * @param local_fields_wanted (out) fields_wanted in local indexes (need to be
 * allocated).
 * @param local_first_deriv (out) Fields (local indexes) corresponding to the
 * first derivative of local_fields_wanted (need to be allocated).
 * @param local_second_deriv (out) Fields (local indexes) corresponding to the
 * second derivative of local_fields_wanted (need to be allocated).
 * @param n_fields_wanted Number of elements in global_fields_wanted,
 * local_fields_wanted and its derivatives.
 * @param type The type of the particle.
 */
void logger_reader_global_to_local(
    const struct logger_reader *reader, const int *global_fields_wanted,
    int *local_fields_wanted, int *local_first_deriv, int *local_second_deriv,
    const int n_fields_wanted, enum part_type type) {

  const struct header *h = &reader->log.header;

  /* Get the correct variables. */
  int n_max = 0;
  int *local_to_global = NULL;
  const char **local_names = NULL;
  switch (type) {
    case swift_type_gas:
      n_max = hydro_logger_field_count;
      local_to_global = hydro_logger_local_to_global;
      local_names = hydro_logger_field_names;
      break;
    case swift_type_dark_matter:
      n_max = gravity_logger_field_count;
      local_to_global = gravity_logger_local_to_global;
      local_names = gravity_logger_field_names;
      break;
    case swift_type_stars:
      n_max = stars_logger_field_count;
      local_to_global = stars_logger_local_to_global;
      local_names = stars_logger_field_names;
      break;
    default:
      error("Particle type not implemented yet.");
  }

  /* Initialize the arrays */
  for (int local = 0; local < n_fields_wanted; local++) {
    local_fields_wanted[local] = -1;
    local_first_deriv[local] = -1;
    local_second_deriv[local] = -1;
  }

  /* Find the corresponding local fields */
  for (int i = 0; i < n_fields_wanted; i++) {
    const int global_field = global_fields_wanted[i];
    const int global_first = h->masks[global_field].reader.first_deriv;
    const int global_second = h->masks[global_field].reader.second_deriv;

    for (int local = 0; local < n_max; local++) {
      /* Check if we have the same field. */
      if (global_field == local_to_global[local]) {
        local_fields_wanted[i] = local;
      }
      /* Check if we have the same first derivative. */
      if (global_first == local_to_global[local]) {
        local_first_deriv[i] = local;
      }
      /* Check if we have the same second derivative. */
      if (global_second == local_to_global[local]) {
        local_second_deriv[i] = local;
      }
    }
  }

  /* Check that we found the fields */
  for (int local = 0; local < n_fields_wanted; local++) {
    if (local_fields_wanted[local] < 0) {
      error("Field %s not found in particle type %s", local_names[local],
            part_type_names[type]);
    }
  }
}

/**
 * @brief Read all the particles from the index file.
 *
 * @param reader The #logger_reader.
 * @param time The requested time for the particle.
 * @param interp_type The type of interpolation.
 * @param global_fields_wanted The fields requested (global index).
 * @param n_fields_wanted Number of field requested.
 * @param output Pointer to the output array. Size: (n_fields_wanted,
 * sum(n_part)).
 * @param n_part Number of particles of each type.
 */
void logger_reader_read_all_particles(struct logger_reader *reader, double time,
                                      enum logger_reader_type interp_type,
                                      const int *global_fields_wanted,
                                      const int n_fields_wanted, void **output,
                                      const uint64_t *n_part) {

  const struct header *h = &reader->log.header;

  /* Allocate temporary memory. */
  /* fields_wanted sorted according to the fields order (local index). */
  int *local_fields_wanted = (int *)malloc(sizeof(int) * n_fields_wanted);
  if (local_fields_wanted == NULL) {
    error("Failed to allocate the array of sorted fields.");
  }

  /* Fields corresponding to the first derivative of fields_wanted (sorted and
   * local index). */
  int *local_first_deriv = malloc(sizeof(int) * n_fields_wanted);
  if (local_first_deriv == NULL) {
    error("Failed to allocate the list of first derivative.");
  }

  /* Fields corresponding to the second derivative of fields_wanted (sorted and
   * local index). */
  int *local_second_deriv = malloc(sizeof(int) * n_fields_wanted);
  if (local_second_deriv == NULL) {
    error("Failed to allocate the list of second derivative.");
  }

  /* Do the hydro. */
  if (n_part[swift_type_gas] != 0) {
    struct index_data *data =
        logger_index_get_data(&reader->index.index, swift_type_gas);

    /* Sort the fields in order to read the correct bits. */
    logger_reader_global_to_local(
        reader, global_fields_wanted, local_fields_wanted, local_first_deriv,
        local_second_deriv, n_fields_wanted, swift_type_gas);

    /* Read the particles */
    for (size_t i = 0; i < n_part[swift_type_gas]; i++) {
      /* Get the offset */
      size_t offset = data[i].offset;

      /* Sort the output into output_single. */
      for (int field = 0; field < n_fields_wanted; field++) {
        const int global = global_fields_wanted[field];
        const int local = local_fields_wanted[field];
        const int first = local_first_deriv[field];
        const int second = local_second_deriv[field];
        void *output_single = output[field] + i * h->masks[global].size;

        /* Read the field. */
        logger_reader_read_field(reader, time, reader->time.time_offset,
                                 interp_type, offset, local, first, second,
                                 output_single, swift_type_gas);
      }
    }
  }

  /* Do the dark matter. */
  if (n_part[swift_type_dark_matter] != 0) {
    struct index_data *data =
        logger_index_get_data(&reader->index.index, swift_type_dark_matter);

    /* Sort the fields in order to read the correct bits. */
    logger_reader_global_to_local(
        reader, global_fields_wanted, local_fields_wanted, local_first_deriv,
        local_second_deriv, n_fields_wanted, swift_type_dark_matter);

    /* Read the particles */
    for (size_t i = 0; i < n_part[swift_type_dark_matter]; i++) {
      /* Get the offset */
      size_t offset = data[i].offset;

      /* Sort the output into output_single. */
      for (int field = 0; field < n_fields_wanted; field++) {
        const int global = global_fields_wanted[field];
        const int local = local_fields_wanted[field];
        const int first = local_first_deriv[field];
        const int second = local_second_deriv[field];
        void *output_single = output[field] + i * h->masks[global].size;

        /* Read the field. */
        logger_reader_read_field(reader, time, reader->time.time_offset,
                                 interp_type, offset, local, first, second,
                                 output_single, swift_type_dark_matter);
      }
    }
  }

  /* Do the stars. */
  if (n_part[swift_type_stars] != 0) {
    struct index_data *data =
        logger_index_get_data(&reader->index.index, swift_type_stars);

    /* Sort the fields in order to read the correct bits. */
    logger_reader_global_to_local(
        reader, global_fields_wanted, local_fields_wanted, local_first_deriv,
        local_second_deriv, n_fields_wanted, swift_type_stars);

    /* Read the particles */
    for (size_t i = 0; i < n_part[swift_type_stars]; i++) {
      /* Get the offset */
      size_t offset = data[i].offset;

      /* Sort the output into output_single. */
      for (int field = 0; field < n_fields_wanted; field++) {
        const int global = global_fields_wanted[field];
        const int local = local_fields_wanted[field];
        const int first = local_first_deriv[field];
        const int second = local_second_deriv[field];
        void *output_single = output[field] + i * h->masks[global].size;

        /* Read the field. */
        logger_reader_read_field(reader, time, reader->time.time_offset,
                                 interp_type, offset, local, first, second,
                                 output_single, swift_type_stars);
      }
    }
  }

  /* Free the memory. */
  free(local_fields_wanted);
  free(local_first_deriv);
  free(local_second_deriv);
}

/**
 * @brief Get the simulation initial time.
 *
 * @param reader The #logger_reader.
 *
 * @return The initial time
 */
double logger_reader_get_time_begin(struct logger_reader *reader) {
  return reader->log.times.records[0].time;
}

/**
 * @brief Get the simulation final time.
 *
 * @param reader The #logger_reader.
 *
 * @return The final time
 */
double logger_reader_get_time_end(struct logger_reader *reader) {
  const size_t ind = reader->log.times.size;
  return reader->log.times.records[ind - 1].time;
}

/**
 * @brief Get the offset of the last timestamp before a given time.
 *
 * @param reader The #logger_reader.
 * @param time The requested time.
 *
 * @return The offset of the timestamp.
 */
size_t logger_reader_get_next_offset_from_time(struct logger_reader *reader,
                                               double time) {
  size_t ind = time_array_get_index_from_time(&reader->log.times, time);
  /* We do not want to have the sentiel */
  if (reader->log.times.size - 2 == ind) {
    ind -= 1;
  }
  return reader->log.times.records[ind + 1].offset;
}

/**
 * @brief Read a record without knowing if it is a particle or a timestamp.
 *
 * WARNING This function asssumes that all the particles are hydro particles.
 * Thus it should be used only for testing the code.
 *
 * @param reader The #logger_reader.
 * @param output The already allocated buffer containing all the fields possible
 * for an hydro particle. (out) The particle if the record is a particle
 * @param time (out) The time if the record is a timestamp.
 * @param is_particle (out) 1 if the record is a particle 0 otherwise.
 * @param offset The offset of the record to read.
 *
 * @return The offset after the record.
 */
size_t logger_reader_read_record(struct logger_reader *reader, void **output,
                                 double *time, int *is_particle,
                                 size_t offset) {

  /* Get a few pointers. */
  const struct header *h = &reader->log.header;
  void *map = reader->log.log.map;

  size_t mask = 0;
  size_t h_offset = 0;

  /* Read the record's mask. */
  map = logger_loader_io_read_mask(h, (char *)map + offset, &mask, &h_offset);

  *is_particle = !(mask & h->timestamp_mask);
  /* The record is a particle. */
  if (*is_particle) {

    size_t offset_tmp = offset;
    for (int i = 0; i < hydro_logger_field_count; i++) {
      offset = logger_particle_read_field(
          reader, offset_tmp, output[i], hydro_logger_local_to_global,
          hydro_logger_field_count, i, &mask, &h_offset);
    }

  }
  /* The record is a timestamp. */
  else {
    integertime_t not_used = 0;
    offset = time_read(&not_used, time, reader, offset);
  }

  return offset;
}
