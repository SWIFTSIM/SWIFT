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
 * @file logger_parameters.h
 * @brief This file contains the functions that deals with the yaml file.
 */

#ifndef LOGGER_LOGGER_PARAMETERS_H
#define LOGGER_LOGGER_PARAMETERS_H

struct logger_reader;

/**
 * @brief Structure containing the simulation's parameters.
 */
struct logger_parameters {

  /* Number of dimension */
  int dimension;

  /* Simulation volume */
  double box_size[3];

  /* Is the volume periodic? */
  int periodic;

  /* Are we doing a cosmological simulation? */
  int with_cosmology;

  /* Number of threads to use in the threadpools */
  int number_threads;
};

void logger_parameters_init(struct logger_parameters *params,
                            const struct logger_reader *reader,
                            int number_threads);

#endif  // LOGGER_LOGGER_PARAMETERS_H
