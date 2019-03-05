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
 * @brief This file contains the C function shown to the external user.
 */
#ifndef __LOGGER_LOGGER_READER_H__
#define __LOGGER_LOGGER_READER_H__

#include "logger_dump.h"
#include "logger_index.h"

/**
 * @brief Main structure of the logger.
 *
 * This structure contains all the variables required for the logger.
 * It should be the only structure that the user see.
 */
struct logger_reader {

  /* Information contained in the index file */
  struct logger_index index;
  
  /* Informations contained in the file header */
  struct logger_dump dump;

  /* Level of verbosity */
  int verbose;
};

void logger_reader_init(struct logger_reader *reader, char *filename, int verbose);
void logger_reader_free(struct logger_reader *reader);

#endif  // __LOGGER_LOGGER_READER_H__
