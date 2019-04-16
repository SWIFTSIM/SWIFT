/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_DEFAULT_BLACK_HOLES_IO_H
#define SWIFT_DEFAULT_BLACK_HOLES_IO_H

#include "black_holes_part.h"
#include "io_properties.h"

/**
 * @brief Specifies which b-particle fields to read from a dataset
 *
 * @param bparts The b-particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
INLINE static void black_holes_read_particles(struct bpart *bparts,
                                              struct io_props *list,
                                              int *num_fields) {

  /* Say how much we want to read */
  *num_fields = 5;

  /* List what we want to read */
  list[0] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                UNIT_CONV_LENGTH, bparts, x);
  list[1] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                UNIT_CONV_SPEED, bparts, v);
  list[2] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                bparts, mass);
  list[3] = io_make_input_field("ParticleIDs", LONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, bparts, id);
  list[4] = io_make_input_field("SmoothingLength", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_LENGTH, bparts, h);
}

/**
 * @brief Specifies which b-particle fields to write to a dataset
 *
 * @param bparts The b-particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
INLINE static void black_holes_write_particles(const struct bpart *bparts,
                                               struct io_props *list,
                                               int *num_fields) {

  /* Say how much we want to write */
  *num_fields = 5;

  /* List what we want to write */
  list[0] = io_make_output_field("Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH,
                                 bparts, x);
  list[1] =
      io_make_output_field("Velocities", FLOAT, 3, UNIT_CONV_SPEED, bparts, v);
  list[2] =
      io_make_output_field("Masses", FLOAT, 1, UNIT_CONV_MASS, bparts, mass);
  list[3] = io_make_output_field("ParticleIDs", LONGLONG, 1, UNIT_CONV_NO_UNITS,
                                 bparts, id);
  list[4] = io_make_output_field("SmoothingLength", FLOAT, 1, UNIT_CONV_LENGTH,
                                 bparts, h);

#ifdef DEBUG_INTERACTIONS_BLACK_HOLES

  list += *num_fields;
  *num_fields += 4;

  list[0] = io_make_output_field("Num_ngb_density", INT, 1, UNIT_CONV_NO_UNITS,
                                 bparts, num_ngb_density);
  list[1] = io_make_output_field("Num_ngb_force", INT, 1, UNIT_CONV_NO_UNITS,
                                 bparts, num_ngb_force);
  list[2] = io_make_output_field("Ids_ngb_density", LONGLONG,
                                 MAX_NUM_OF_NEIGHBOURS_BLACK_HOLES,
                                 UNIT_CONV_NO_UNITS, bparts, ids_ngbs_density);
  list[3] = io_make_output_field("Ids_ngb_force", LONGLONG,
                                 MAX_NUM_OF_NEIGHBOURS_BLACK_HOLES,
                                 UNIT_CONV_NO_UNITS, bparts, ids_ngbs_force);
#endif
}

#endif /* SWIFT_DEFAULT_BLACK_HOLES_IO_H */
