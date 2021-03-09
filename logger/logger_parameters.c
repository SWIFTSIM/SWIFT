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
#include "logger_parameters.h"

/* Include the logger */
#include "logger_reader.h"
#include "logger_tools.h"

/* Include swift */
#include "parser.h"

/**
 * @brief Initialize the parameters from the yaml file.
 *
 * @param reader The #logger_reader.
 * @param basename The basename of the logger files.
 * @param verbose The verbose level.
 */
void logger_parameters_init(struct logger_parameters *params,
                            const struct logger_reader *reader) {

  /* Generate filename */
  char filename[STRING_SIZE];
  strcpy(filename, reader->basename);
  size_t len = strlen(filename);
  /* Remove the rank number */
  filename[len - 5] = '\0';
  strcat(filename, ".yml");

  /* Initialize the parser */
  struct swift_params swift_params;
  parser_read_file(filename, &swift_params);

#ifdef SWIFT_DEBUG_CHECKS
  /* Print the parameters */
  if (reader->verbose > 0) parser_print_params(&swift_params);
#endif

  /* Read the number of dimension */
  params->dimension = parser_get_param_int(&swift_params, "Header:Dimension");

  /* Read the box size */
  parser_get_param_double_array(&swift_params, "Header:BoxSize", 3,
                                params->box_size);

  /* Read the periodicity */
  params->periodic = parser_get_param_int(&swift_params, "Header:Periodic");

  /* Read if we are running with cosmology */
  params->with_cosmology =
      parser_get_param_int(&swift_params, "Policy:cosmological integration");
}
