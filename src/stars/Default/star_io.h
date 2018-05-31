/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_DEFAULT_STAR_IO_H
#define SWIFT_DEFAULT_STAR_IO_H

#include "io_properties.h"

/**
 * @brief Specifies which s-particle fields to read from a dataset
 *
 * @param sparts The s-particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
INLINE static void star_read_particles(struct spart* sparts,
                                       struct io_props* list, int* num_fields) {

  /* Say how much we want to read */
  *num_fields = 4;

  /* List what we want to read */
  list[0] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                UNIT_CONV_LENGTH, sparts, x);
  list[1] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                UNIT_CONV_SPEED, sparts, v);
  list[2] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                sparts, mass);
  list[3] = io_make_input_field("ParticleIDs", LONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, sparts, id);
}

/**
 * @brief Specifies which s-particle fields to write to a dataset
 *
 * @param sparts The s-particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
INLINE static void star_write_particles(const struct spart* sparts,
                                        struct io_props* list,
                                        int* num_fields) {

  /* Say how much we want to read */
  *num_fields = 4;

  /* List what we want to read */
  list[0] = io_make_output_field("Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH,
                                 sparts, x);
  list[1] =
      io_make_output_field("Velocities", FLOAT, 3, UNIT_CONV_SPEED, sparts, v);
  list[2] =
      io_make_output_field("Masses", FLOAT, 1, UNIT_CONV_MASS, sparts, mass);
  list[3] = io_make_output_field("ParticleIDs", LONGLONG, 1, UNIT_CONV_NO_UNITS,
                                 sparts, id);
}

#endif /* SWIFT_DEFAULT_STAR_IO_H */
