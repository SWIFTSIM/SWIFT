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
#ifndef SWIFT_DEFAULT_GRAVITY_IO_H
#define SWIFT_DEFAULT_GRAVITY_IO_H

#include "io_properties.h"

void convert_gpart_pos(const struct engine* e, const struct gpart* gp,
                       double* ret) {

  if (e->s->periodic) {
    ret[0] = box_wrap(gp->x[0], 0.0, e->s->dim[0]);
    ret[1] = box_wrap(gp->x[1], 0.0, e->s->dim[1]);
    ret[2] = box_wrap(gp->x[2], 0.0, e->s->dim[2]);
  } else {
    ret[0] = gp->x[0];
    ret[1] = gp->x[1];
    ret[2] = gp->x[2];
  }
}

/**
 * @brief Specifies which g-particle fields to read from a dataset
 *
 * @param gparts The g-particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
void darkmatter_read_particles(struct gpart* gparts, struct io_props* list,
                               int* num_fields) {

  /* Say how much we want to read */
  *num_fields = 4;

  /* List what we want to read */
  list[0] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                UNIT_CONV_LENGTH, gparts, x);
  list[1] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                UNIT_CONV_SPEED, gparts, v_full);
  list[2] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                gparts, mass);
  list[3] = io_make_input_field("ParticleIDs", ULONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, gparts, id_or_neg_offset);
}

/**
 * @brief Specifies which g-particle fields to write to a dataset
 *
 * @param gparts The g-particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
void darkmatter_write_particles(const struct gpart* gparts,
                                struct io_props* list, int* num_fields) {

  /* Say how much we want to read */
  *num_fields = 5;

  /* List what we want to read */
  list[0] = io_make_output_field_convert_gpart(
      "Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH, gparts, convert_gpart_pos);
  list[1] = io_make_output_field("Velocities", FLOAT, 3, UNIT_CONV_SPEED,
                                 gparts, v_full);
  list[2] =
      io_make_output_field("Masses", FLOAT, 1, UNIT_CONV_MASS, gparts, mass);
  list[3] = io_make_output_field("ParticleIDs", ULONGLONG, 1,
                                 UNIT_CONV_NO_UNITS, gparts, id_or_neg_offset);
  list[4] = io_make_output_field("Acceleration", FLOAT, 3,
                                 UNIT_CONV_ACCELERATION, gparts, a_grav);
}

#endif /* SWIFT_DEFAULT_GRAVITY_IO_H */
