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
 * @file logger_reader.h
 * @brief This file contains the C functions shown to the external user.
 *
 * Here is a quick summary of our different elements:
 *
 * The logger is a time adaptive way to write snapshots.
 * It consists of a set of files: the log file, the parameter file and the index
 * files.
 *
 * The <b>parameter file</b> contains all the information related to the code
 * (e.g. boxsize, scheme, ...).
 *
 * The <b>index files</b> are not mandatory files that indicates the position of
 * the particles in the log file at a given time step. They are useful to
 * speedup the reading.
 *
 * The <b>log file</b> consists in a large file where the particles are logged
 * one after the other. It contains a <b>log file header</b> at the beginning of
 * the file and a large collection of <b>records</b>.
 *
 * The records are logged one after the other and each contains a <b>record
 * header</b> and then a list of <b>named entries</b>. In the record header, a
 * <b>mask</b> is provided that corresponds to the type of named entries present
 * in this record. It also contains the <b>offset</b> to the previous or next
 * record for this particle.
 */

#ifndef LOGGER_LOGGER_READER_H
#define LOGGER_LOGGER_READER_H

#include "logger_index.h"
#include "logger_loader_io.h"
#include "logger_logfile.h"
#include "logger_parameters.h"
#include "logger_particle.h"

/**
 * @brief Main structure of the logger.
 *
 * This structure contains all the variables required for the logger.
 * It should be the only structure that the user see.
 *
 * It is initialized with #logger_reader_init and freed with
 * #logger_reader_free.
 */
struct logger_reader {

  /* Base name of the files */
  char basename[STRING_SIZE];

  struct {
    /* Information contained in the previous index file. */
    struct logger_index index_prev;

    /* Information contained in the next index file. */
    struct logger_index index_next;

    /* Number of index files */
    int n_files;

    /* Time of each index file */
    double *times;

    /* Integer time of each index file */
    integertime_t *int_times;
  } index;

  /* Informations contained in the file header. */
  struct logger_logfile log;

  /* Information about the current time */
  struct {
    /* Double time */
    double time;

    /* Integer time */
    integertime_t int_time;

    /* Offset of the chunk */
    size_t time_offset;

    /* Index of the element in the time array */
    size_t index;
  } time;

  /* Information from the yaml file */
  struct logger_parameters params;

  /* Level of verbosity. */
  int verbose;
};

enum logger_reader_event {
  logger_reader_event_null,    /* No event */
  logger_reader_event_deleted, /* Particle has been deleted */
  logger_reader_event_stars, /* The particle has been transformed into a star */
};

void logger_reader_init_index(struct logger_reader *reader);
void logger_reader_init(struct logger_reader *reader, const char *basename,
                        int verbose, int number_threads);
void logger_reader_free(struct logger_reader *reader);

void logger_reader_set_time(struct logger_reader *reader, double time);

double logger_reader_get_time_begin(struct logger_reader *reader);
double logger_reader_get_time_end(struct logger_reader *reader);
size_t logger_reader_get_next_offset_from_time(struct logger_reader *reader,
                                               double time);
void logger_reader_get_number_particles(struct logger_reader *reader,
                                        uint64_t *n_parts,
                                        const int *read_types);

void logger_reader_read_all_particles(struct logger_reader *reader, double time,
                                      enum logger_reader_type interp_type,
                                      const int *id_masks_wanted,
                                      const int n_mask_wanted, void **output,
                                      const uint64_t *n_part);
void logger_reader_read_particles_from_ids(
    struct logger_reader *reader, double time,
    enum logger_reader_type interp_type, const int *id_masks_wanted,
    const int n_mask_wanted, void **output, uint64_t *n_part, long long **ids);
size_t logger_reader_read_record(struct logger_reader *reader, void **output,
                                 double *time, int *is_particle, size_t offset);

#endif  // LOGGER_LOGGER_READER_H
