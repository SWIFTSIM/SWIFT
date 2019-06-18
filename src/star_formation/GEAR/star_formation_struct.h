/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Folkert Nobels (nobels@strw.leidenuniv.nl)
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
#ifndef SWIFT_GEAR_STAR_FORMATION_STRUCT_H
#define SWIFT_GEAR_STAR_FORMATION_STRUCT_H

/**
 * @brief Star-formation-related properties stored in the extended particle
 * data.
 */
struct star_formation_xpart_data {};


struct star_formation_part_data {
  // TODO move it to the pressure floor
  /*! Estimation of local turbulence (squared) */
  float sigma2;
};

/**
 * @brief Global star formation properties
 */
struct star_formation {

  // TODO move it to pressure floor
  /*! Number of particle required to resolved the Jeans criterion (at power 2/3) */
  float n_jeans_2_3;

  /*! Maximal temperature for forming a star */
  float maximal_temperature;

  /*! Star formation efficiency */
  float star_formation_efficiency;
};

#endif /* SWIFT_GEAR_STAR_FORMATION_STRUCT_H */
