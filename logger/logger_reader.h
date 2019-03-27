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
 * It consists of a set of files: the log file, the parameter file and the index files.
 *
 * The <b>parameter file</b> contains all the information related to the code (e.g. boxsize).
 *
 * The <b>index files</b> are not mandatory files that indicates the position of the particles in
 * the log file at a given time step. They are useful to speedup the reading.
 *
 * The <b>log file</b> consists in a large file where the particles are logged one after the other.
 * It contains a <b>log file header</b> at the beginning of the file and a large collection of <b>records</b>.
 *
 * The records are logged one after the other and each contains a <b>record header</b> and then a list of <b>named entries</b>.
 * In the record header, a <b>mask</b> is provided that corresponds to the type of named entries present in this record.
 * It also contains the <b>offset</b> to the previous or next record for this particle.
 */

#ifndef __LOGGER_LOGGER_READER_H__
#define __LOGGER_LOGGER_READER_H__

#include "logger_logfile.h"
#include "logger_index.h"

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

  /* Information contained in the index file */
  struct logger_index index;
  
  /* Informations contained in the file header */
  struct logger_logfile log;

  /* Level of verbosity */
  int verbose;
};

void logger_reader_init(struct logger_reader *reader, char *filename, int verbose);
void logger_reader_free(struct logger_reader *reader);

#endif  // __LOGGER_LOGGER_READER_H__
