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

#include "logger_reader.h"

/**
 * @brief Initialize the reader.
 *
 * @param reader The #logger_reader.
 * @param filename The log filename.
 * @param verbose The verbose level.
 */
void logger_reader_init(struct logger_reader *reader, char *filename,
                        int verbose) {
  if (verbose > 1) message("Initializing the reader.");

  /* Initialize the reader variables. */
  reader->verbose = verbose;

  /* Initialize the log file. */
  logger_logfile_init_from_file(&reader->log, filename, reader,
                                /* only_header */ 0);

  if (verbose > 1) message("Initialization done.");
}

/**
 * @brief Free the reader.
 *
 * @param reader The #logger_reader.
 */
void logger_reader_free(struct logger_reader *reader) {
  /* Free the log. */
  logger_logfile_free(&reader->log);
}

/**
 * @brief Read a record (timestamp or particle)
 *
 * @param reader The #logger_reader.
 * @param lp (out) The #logger_particle (if the record is a particle).
 * @param time (out) The time read (if the record is a timestamp).
 * @param is_particle Is the record a particle (or a timestamp)?
 * @param offset The offset in the file.
 *
 * @return The offset after this record.
 */
size_t reader_read_record(struct logger_reader *reader,
                          struct logger_particle *lp, double *time,
                          int *is_particle, size_t offset) {

  struct logger_logfile *log = &reader->log;

  /* Read mask to find out if timestamp or particle. */
  size_t mask = 0;
  logger_loader_io_read_mask(&log->header, log->log.map + offset, &mask, NULL);

  /* Check if timestamp or not. */
  int ind = header_get_field_index(&log->header, "timestamp");
  if (ind == -1) {
    error("File header does not contain a mask for time.");
  }
  if (log->header.masks[ind].mask == mask) {
    *is_particle = 0;
    integertime_t int_time = 0;
    offset = time_read(&int_time, time, reader, offset);
  } else {
    *is_particle = 1;
    offset =
        logger_particle_read(lp, reader, offset, *time, logger_reader_const);
  }

  return offset;
}
