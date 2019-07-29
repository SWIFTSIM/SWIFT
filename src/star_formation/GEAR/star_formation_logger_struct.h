/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
 *               2019 Fabien Jeanquartier (fabien.jeanquartier@epfl.ch)
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
#ifndef SWIFT_GEAR_STAR_FORMATION_LOGGER_STRUCT_H
#define SWIFT_GEAR_STAR_FORMATION_LOGGER_STRUCT_H

/**
 * Structure containing the star formation information from the cells.
 */
struct star_formation_history {
  /*! Stellar mass created in the current timestep */
  float new_stellar_mass;

  /*! Number of stars created in this timestep */
  long int number_new_stars;
};

/**
 * Structure containing the global star formation information (including time
 * integrated variables).
 */
struct star_formation_history_accumulator {
  /*! Total stellar mass from the begining of the simulation */
  float total_stellar_mass;

  /*! Total number of stars */
  long int total_number_stars;

  /*! Stellar mass created in the current timestep */
  float new_stellar_mass;

  /*! Number of stars created in this timestep */
  long int number_new_stars;
};

#endif /* SWIFT_GEAR_STAR_FORMATION_LOGGER_STRUCT_H */
