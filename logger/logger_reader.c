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
  logger_index_init(&reader->index.index_prev, reader);
  logger_index_init(&reader->index.index_next, reader);

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
    error_python("Failed to allocate the list of times");
  }
  if ((reader->index.int_times =
           (integertime_t *)malloc(count * sizeof(integertime_t))) == NULL) {
    error_python("Failed to allocate the list of times");
  }

  /* Get the information contained in the headers */
  for (int i = 0; i < reader->index.n_files; i++) {
    char filename[STRING_SIZE + 50];
    sprintf(filename, "%s_%04i.index", reader->basename, i);

    /* Read the header */
    logger_index_read_header(&reader->index.index_prev, filename);

    /* Save the required information */
    reader->index.times[i] = reader->index.index_prev.time;
    reader->index.int_times[i] = reader->index.index_prev.integer_time;
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
    logger_index_free(&reader->index.index_prev);
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
  char filename_prev[STRING_SIZE + 50];
  sprintf(filename_prev, "%s_%04u.index", reader->basename, left);
  char filename_next[STRING_SIZE + 50];
  sprintf(filename_next, "%s_%04u.index", reader->basename, left + 1);

  /* Check if the file is already mapped */
  if (reader->index.index_prev.index.map != NULL) {
    logger_index_free(&reader->index.index_prev);
  }
  if (reader->index.index_next.index.map != NULL) {
    logger_index_free(&reader->index.index_next);
  }

  /* Read the file */
  logger_index_read_header(&reader->index.index_prev, filename_prev);
  logger_index_map_file(&reader->index.index_prev, filename_prev,
                        /* sorted */ 1);

  logger_index_read_header(&reader->index.index_next, filename_next);
  logger_index_map_file(&reader->index.index_next, filename_next,
                        /* sorted */ 1);

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
 * @brief Count the number of new particles at the time set since last index
 * file.
 *
 * @param reader The #logger_reader.
 * @param part_type The index corresponding to the particle type.
 *
 * @param The number of new particles.
 */
size_t logger_reader_count_number_new_particles(struct logger_reader *reader,
                                                enum part_type part_type) {

  const size_t threshold = reader->time.time_offset;

  /* Get the history of created particles. */
  struct index_data *data =
      logger_index_get_created_history(&reader->index.index_next, part_type);

  /* Do we have any new particle? */
  if (reader->index.index_next.nparts_created[part_type] == 0 ||
      threshold < data[0].offset) {
    return 0;
  }

  /* Perform a binary search */
  size_t right = reader->index.index_next.nparts_created[part_type] - 1;
  size_t left = 0;
  while (left <= right) {
    size_t m = (left + right) / 2;
    if (data[m].offset < threshold) {
      left = m + 1;
    } else if (data[m].offset > threshold) {
      right = m - 1;
    } else {
      return m;
    }
  }

  /* Compute the return value. */
  size_t ret = (left + right) / 2;

  /* Check if the binary search is correct. */
  if (data[ret].offset > threshold) {
    error_python("The binary search has failed.");
  }

  /* We need the count, not the index. */
  return ret + 1;
}

/**
 * @brief Count the number of removed particles at the time set since last index
 * file.
 *
 * @param reader The #logger_reader.
 * @param n_type (output) The number of particle type possible.
 * @param part_type The index corresponding to the particle type.
 *
 * @param The number of removed particles.
 */
size_t logger_reader_count_number_removed_particles(
    struct logger_reader *reader, enum part_type part_type) {

  const size_t threshold = reader->time.time_offset;

  /* Get the history of created particles. */
  struct index_data *data =
      logger_index_get_removed_history(&reader->index.index_next, part_type);

  /* Do we have any new particle? */
  if (reader->index.index_next.nparts_removed[part_type] == 0 ||
      threshold < data[0].offset) {
    return 0;
  }

  /* Perform a binary search */
  size_t right = reader->index.index_next.nparts_removed[part_type] - 1;
  size_t left = 0;
  while (left <= right) {
    size_t m = (left + right) / 2;
    if (data[m].offset < threshold) {
      left = m + 1;
    } else if (data[m].offset > threshold) {
      right = m - 1;
    } else {
      return m;
    }
  }

  /* Compute the return value. */
  size_t ret = (left + right) / 2;

  /* Check if the binary search is correct. */
  if (data[ret].offset > threshold) {
    error_python("The binary search has failed.");
  }

  /* We need the count, not the index. */
  return ret + 1;
}

/**
 * @brief Provides the number of particle (per type) from the index file.
 *
 * @param reader The #logger_reader.
 * @param n_parts (out) Number of particles at the time set in the reader.
 * @param read_types Should we read this type of particles?
 */
void logger_reader_get_number_particles(struct logger_reader *reader,
                                        uint64_t *n_parts,
                                        const int *read_types) {
  for (enum part_type i = (enum part_type)0; i < swift_type_count; i++) {
    /* Should we skip this type of particles? */
    if (read_types[i] == 0) {
      n_parts[i] = 0;
      continue;
    }

    /* Count the number of particles present in the last index file. */
    const uint64_t count_prev = reader->index.index_prev.nparts[i];
    /* Count the number of particles that have been created since last index. */
    const uint64_t count_new =
        logger_reader_count_number_new_particles(reader, i);
    /* Count the number of particles that have been removed since last index. */
    const uint64_t count_removed =
        logger_reader_count_number_removed_particles(reader, i);
    n_parts[i] = (count_prev + count_new) - count_removed;
  }
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
 * @param field_wanted The field to read.
 * @param all_fields The list of all the fields in the particle type.
 * @param all_fields_count The number of elements in all_fields.
 * @param output Pointers to the output array
 * @param type The #part_type.
 *
 * @return Is the particle removed from the logfile?
 */
int logger_reader_read_field(struct logger_reader *reader, double time,
                             size_t offset_time,
                             enum logger_reader_type interp_type,
                             const size_t offset_last_full_record,
                             const struct field_information *field_wanted,
                             const struct field_information *all_fields,
                             const int all_fields_count, void *output,
                             enum part_type type) {

  const struct header *h = &reader->log.header;
  size_t offset = offset_last_full_record;

  /* Get the indexes */
  const int global = field_wanted->global_index;
  const int global_first = field_wanted->global_index_first;
  const int global_second = field_wanted->global_index_second;

  /* Check if the offsets are correct. */
  if (offset > offset_time) {
    error_python("Last offset is larger than the requested time.");
  }

  /* Get the masks. */
  struct mask_data *mask_field = &h->masks[global];
  struct mask_data *mask_first = NULL;
  if (global_first >= 0) {
    mask_first = &h->masks[global_first];
  }
  struct mask_data *mask_second = NULL;
  if (global_second >= 0) {
    mask_second = &h->masks[global_second];
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

    /* Is the particle removed from the logfile? */
    if (mask & h->masks[logger_index_special_flags].mask) {
      int data = 0;
      enum logger_special_flags flag = logger_particle_read_special_flag(
          reader, offset, &mask, &h_offset, &data);
      if (flag == logger_flag_change_type || flag == logger_flag_mpi_exit ||
          flag == logger_flag_delete) {
        return 1;
      }
    }

    /* Check if there is a next record (avoid infinite loop). */
    if (h_offset == 0) {
      error_python("No next record found.");
    }

    /* Check if the field is present */
    if (mask & mask_field->mask) {
      offset_before = offset;
    }

    /* Go to the next record. */
    offset += h_offset;
  }

  /* Read the field */
  logger_particle_read_field(reader, offset_before, output, all_fields,
                             all_fields_count, global, &mask, &h_offset);

  /* Deal with the first derivative. */
  const size_t size_first = mask_first == NULL ? 0 : mask_first->size;
  char first_deriv[size_first];

  int first_found = mask_first != NULL && mask & mask_first->mask;
  if (first_found) {
    /* Read the first derivative */
    logger_particle_read_field(reader, offset_before, first_deriv, all_fields,
                               all_fields_count, global_first, &mask,
                               &h_offset);
  }

  /* Deal with the second derivative. */
  const size_t size_second = mask_second == NULL ? 0 : mask_second->size;
  char second_deriv[size_second];

  int second_found = mask_second != NULL && mask & mask_second->mask;
  if (second_found) {
    /* Read the first derivative */
    logger_particle_read_field(reader, offset_before, second_deriv, all_fields,
                               all_fields_count, global_second, &mask,
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
        error_python("There is no record after the current one");
      }

      /* Go to the next record. */
      offset += h_offset;
    }

    /* Output after the requested time. */
    char output_after[mask_field->size];

    /* Read the field */
    logger_particle_read_field(reader, offset, output_after, all_fields,
                               all_fields_count, global, &mask, &h_offset);

    /* Deal with the first derivative. */
    char first_deriv_after[size_first];

    /* Did we find the derivative before and in this record? */
    first_found = mask_first != NULL && first_found && mask & mask_first->mask;
    if (first_found) {
      /* Read the first derivative */
      logger_particle_read_field(reader, offset, first_deriv_after, all_fields,
                                 all_fields_count, global_first, &mask,
                                 &h_offset);
    }

    /* Deal with the second derivative. */
    char second_deriv_after[size_second];

    /* Did we find the derivative before and in this record? */
    second_found =
        mask_second != NULL && second_found && mask & mask_second->mask;
    if (second_found) {
      /* Read the second derivative */
      logger_particle_read_field(reader, offset, second_deriv_after, all_fields,
                                 all_fields_count, global_second, &mask,
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
    logger_particle_interpolate_field(time_before, &before, time_after, &after,
                                      output, time, field_wanted, type);
  }
  return 0;
}

/**
 * @brief Convert the fields from global indexes to local.
 *
 * @param reader The #logger_reader.
 * @param global_fields_wanted The fields to convert.
 * @param fields_wanted (out) Fields wanted (need to be allocated).
 * @param n_fields_wanted Number of elements in global_fields_wanted,
 * local_fields_wanted and its derivatives.
 * @param all_fields The list of all the fields.
 * @param all_fields_count The number of elements in all_fields.
 * @param type The type of the particle.
 */
void logger_reader_get_fields_wanted(const struct logger_reader *reader,
                                     const int *global_fields_wanted,
                                     struct field_information *fields_wanted,
                                     const int n_fields_wanted,
                                     struct field_information *all_fields,
                                     int all_fields_count,
                                     enum part_type type) {

  const struct header *h = &reader->log.header;

  /* Mark fields_wanted in order to check that the fields are found. */
  for (int i = 0; i < n_fields_wanted; i++) {
    fields_wanted[i].local_index = -1;
  }

  /* Find the corresponding fields */
  for (int i = 0; i < n_fields_wanted; i++) {
    const int global = global_fields_wanted[i];

    for (int j = 0; j < all_fields_count; j++) {
      /* Copy the structure if the field is found.  */
      if (all_fields[j].global_index == global) {
        fields_wanted[i] = all_fields[j];
        break;
      }
    }
  }

  /* Check that we found the fields */
  for (int i = 0; i < n_fields_wanted; i++) {
    if (fields_wanted[i].local_index < 0) {
      const int global = global_fields_wanted[i];
      const char *name = h->masks[global].name;
      error_python("Field %s not found in particle type %s", name,
                   part_type_names[type]);
    }
  }
}

/**
 * @brief Read all the particles of a given type from the index file.
 *
 * @param reader The #logger_reader.
 * @param time The requested time for the particle.
 * @param interp_type The type of interpolation.
 * @param global_fields_wanted The fields requested (global index).
 * @param n_fields_wanted Number of field requested.
 * @param output Pointer to the output array. Size: (n_fields_wanted,
 * sum(n_part)).
 * @param n_part Number of particles of each type.
 * @param type The particle type
 */
void logger_reader_read_all_particles_single_type(
    struct logger_reader *reader, double time,
    enum logger_reader_type interp_type, const int *global_fields_wanted,
    const int n_fields_wanted, void **output, const uint64_t *n_part,
    enum part_type type) {

  const struct header *h = &reader->log.header;

  /* Count the number of previous parts for the shift in output */
  uint64_t prev_npart = 0;
  for (int i = 0; i < type; i++) {
    prev_npart += n_part[i];
  }

  /* Allocate temporary memory. */
  struct field_information *fields_wanted = (struct field_information *)malloc(
      sizeof(struct field_information) * n_fields_wanted);
  if (fields_wanted == NULL) {
    error_python("Failed to allocate the field information.");
  }

  struct index_data *data =
      logger_index_get_data(&reader->index.index_prev, type);
  struct index_data *data_created =
      logger_index_get_created_history(&reader->index.index_next, type);

  /* Get the list of fields. */
  const int all_fields_count = tools_get_number_fields(type);
  struct field_information *all_fields = (struct field_information *)malloc(
      all_fields_count * sizeof(struct field_information));
  tools_get_list_fields(all_fields, type, h);

  /* Convert fields into the local array. */
  logger_reader_get_fields_wanted(reader, global_fields_wanted, fields_wanted,
                                  n_fields_wanted, all_fields, all_fields_count,
                                  type);

  size_t current_in_index = 0;
  int reading_history = 0;
  const size_t size_index = reader->index.index_prev.nparts[type];
  const size_t size_history =
      logger_reader_count_number_new_particles(reader, type);

  /* Read the particles */
  for (size_t i = 0; i < n_part[type]; i++) {
    int particle_removed = 1;
    /* Do it until finding a particle not removed. */
    while (particle_removed) {
      /* Should we start to read the history? */
      if (!reading_history && current_in_index == size_index) {
        current_in_index = 0;
        reading_history = 1;
      }

      /* Check if we still have some particles available. */
      if (reading_history && current_in_index == size_history) {
        error_python("The logger was not able to find enough particles.");
      }

      /* Get the offset */
      size_t offset = reading_history ? data_created[current_in_index].offset
                                      : data[current_in_index].offset;

      /* Loop over each field. */
      for (int field = 0; field < n_fields_wanted; field++) {
        const int global = fields_wanted[field].global_index;
        void *output_single =
            (char *)output[field] + (i + prev_npart) * h->masks[global].size;

        /* Read the field. */
        particle_removed = logger_reader_read_field(
            reader, time, reader->time.time_offset, interp_type, offset,
            &fields_wanted[field], all_fields, all_fields_count, output_single,
            type);

        /* Should we continue to read the fields of this particle? */
        if (particle_removed) {
          break;
        }
      }
      current_in_index++;
    }
  }

  /* Free the memory */
  free(all_fields);
  free(fields_wanted);
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

  /* Read the gas */
  if (n_part[swift_type_gas] != 0) {
    logger_reader_read_all_particles_single_type(
        reader, time, interp_type, global_fields_wanted, n_fields_wanted,
        output, n_part, swift_type_gas);
  }

  /* Read the dark matter. */
  if (n_part[swift_type_dark_matter] != 0) {
    logger_reader_read_all_particles_single_type(
        reader, time, interp_type, global_fields_wanted, n_fields_wanted,
        output, n_part, swift_type_dark_matter);
  }
  /* Read the dark matter background. */
  if (n_part[swift_type_dark_matter_background] != 0) {
    logger_reader_read_all_particles_single_type(
        reader, time, interp_type, global_fields_wanted, n_fields_wanted,
        output, n_part, swift_type_dark_matter_background);
  }
  /* Read the stars. */
  if (n_part[swift_type_stars] != 0) {
    logger_reader_read_all_particles_single_type(
        reader, time, interp_type, global_fields_wanted, n_fields_wanted,
        output, n_part, swift_type_stars);
  }
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
  char *map = reader->log.log.map;

  size_t mask = 0;
  size_t h_offset = 0;

  /* Read the record's mask. */
  map = logger_loader_io_read_mask(h, map + offset, &mask, &h_offset);

  *is_particle = !(mask & h->timestamp_mask);
  /* The record is a particle. */
  if (*is_particle) {
    /* Get the list of fields. */
    const int all_fields_count = tools_get_number_fields(swift_type_gas);
    struct field_information *all_fields = (struct field_information *)malloc(
        all_fields_count * sizeof(struct field_information));
    tools_get_list_fields(all_fields, swift_type_gas, h);

    size_t offset_tmp = offset;
    for (int i = 0; i < all_fields_count; i++) {
      offset = logger_particle_read_field(
          reader, offset_tmp, output[i], all_fields, all_fields_count,
          all_fields[i].global_index, &mask, &h_offset);
    }

  }
  /* The record is a timestamp. */
  else {
    integertime_t not_used = 0;
    offset = time_read(&not_used, time, reader, offset);
  }

  return offset;
}
