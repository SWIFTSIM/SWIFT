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
void logger_reader_init(struct logger_reader *reader, char *filename, int verbose) {
  if (verbose > 1)
    message("Initializing the reader");
  /* Initialize the reader variables */
  reader->verbose = verbose;

  /* Initialize the log */
  logger_logfile_init(&reader->log, filename, reader);

  if (verbose > 1)
    message("Initialization done.");
}

/**
 * @brief Free the reader.
 *
 * @param reader The #logger_reader.
 */
void logger_reader_free(struct logger_reader *reader) {
  /* Free the log */
  logger_logfile_free(&reader->log);
}
