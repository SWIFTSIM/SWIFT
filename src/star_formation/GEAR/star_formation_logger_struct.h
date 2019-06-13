/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_EAGLE_STAR_FORMATION_LOGGER_STRUCT_H
#define SWIFT_EAGLE_STAR_FORMATION_LOGGER_STRUCT_H

/* Starformation history struct */
struct star_formation_history {
  /*! Total current star formation rate (new mass divided by timestep) */
  float new_sfr;
  /*! Total stellar mass (no reset) */
  float total_stellar_mass;
  /*! Stellar mass created in the current timestep */
  float stellar_mass;
  /*! Number of stars created in this timestep */
  long int number_of_stars;
  /*! Total number of stars */
  long int total_number_of_stars;
};

#endif /* SWIFT_EAGLE_STAR_FORMATION_LOGGER_STRUCT_H */
