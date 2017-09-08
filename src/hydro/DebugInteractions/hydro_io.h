/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2017 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_DEBUG_INTERACTIONS_HYDRO_IO_H
#define SWIFT_DEBUG_INTERACTIONS_HYDRO_IO_H

/**
 * @file DebugInteractions/hydro_io.h
 * @brief Empty SPH implementation used solely to test the SELF/PAIR routines.
 */

#include "adiabatic_index.h"
#include "hydro.h"
#include "io_properties.h"
#include "kernel_hydro.h"

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
void hydro_read_particles(struct part* parts, struct io_props* list,
                          int* num_fields) {

  *num_fields = 4;

  /* List what we want to read */
  list[0] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                UNIT_CONV_LENGTH, parts, x);
  list[1] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                UNIT_CONV_SPEED, parts, v);
  list[2] = io_make_input_field("SmoothingLength", FLOAT, 1, COMPULSORY,
                                UNIT_CONV_LENGTH, parts, h);
  list[3] = io_make_input_field("ParticleIDs", ULONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, parts, id);
}

float convert_u(struct engine* e, struct part* p) {

  return hydro_get_internal_energy(p);
}

float convert_P(struct engine* e, struct part* p) {

  return hydro_get_pressure(p);
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
void hydro_write_particles(struct part* parts, struct io_props* list,
                           int* num_fields) {

  *num_fields = 8;

  /* List what we want to write */
  list[0] = io_make_output_field("Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH,
                                 parts, x);
  list[1] =
      io_make_output_field("Velocities", FLOAT, 3, UNIT_CONV_SPEED, parts, v);
  list[2] = io_make_output_field("SmoothingLength", FLOAT, 1, UNIT_CONV_LENGTH,
                                 parts, h);
  list[3] = io_make_output_field("ParticleIDs", ULONGLONG, 1,
                                 UNIT_CONV_NO_UNITS, parts, id);
  list[4] = io_make_output_field("Num_ngb_density", INT, 1, UNIT_CONV_NO_UNITS,
                                 parts, num_ngb_density);
  list[5] = io_make_output_field("Num_ngb_force", INT, 1, UNIT_CONV_NO_UNITS,
                                 parts, num_ngb_force);
  list[6] = io_make_output_field("Ids_ngb_density", LONGLONG, 256,
                                 UNIT_CONV_NO_UNITS, parts, ids_ngbs_density);
  list[7] = io_make_output_field("Ids_ngb_force", LONGLONG, 256,
                                 UNIT_CONV_NO_UNITS, parts, ids_ngbs_force);
}

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grpsph The HDF5 group in which to write
 */
void writeSPHflavour(hid_t h_grpsph) {

  /* Viscosity and thermal conduction */
  io_write_attribute_s(h_grpsph, "Thermal Conductivity Model", "None");
  io_write_attribute_s(h_grpsph, "Viscosity Model", "None");
}

/**
 * @brief Are we writing entropy in the internal energy field ?
 *
 * @return 1 if entropy is in 'internal energy', 0 otherwise.
 */
int writeEntropyFlag() { return 0; }

#endif /* SWIFT_DEBUG_INTERACTIONS_HYDRO_IO_H */
