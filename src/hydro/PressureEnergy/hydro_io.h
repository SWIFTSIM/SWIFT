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
#ifndef SWIFT_PRESSURE_ENERGY_HYDRO_IO_H
#define SWIFT_PRESSURE_ENERGY_HYDRO_IO_H

/**
 * @file PressureEnergy/hydro_io.h
 * @brief Pressure-Energy implementation of SPH (i/o routines)
 *
 * The thermal variable is the energy (u) and the pressure is smoothed over
 * contact discontinuities to prevent spurious surface tension.
 *
 * Follows equations (17) and (18) of Hopkins, P., MNRAS, 2013,
 * Volume 428, Issue 4, pp. 2840-2856 with a simple Balsara viscosity term.
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

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
void hydro_write_particles(struct part* parts, struct io_props* list,
                           int* num_fields) {

  *num_fields = 4;

  /* List what we want to write */
  list[0] = io_make_output_field("Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH,
                                 parts, x);
  list[1] =
      io_make_output_field("Velocities", FLOAT, 3, UNIT_CONV_SPEED, parts, v);
  list[2] = io_make_output_field("SmoothingLength", FLOAT, 1, UNIT_CONV_LENGTH,
                                 parts, h);
  list[3] = io_make_output_field("ParticleIDs", ULONGLONG, 1,
                                 UNIT_CONV_NO_UNITS, parts, id);
}

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grpsph The HDF5 group in which to write
 */
void writeSPHflavour(hid_t h_grpsph) {}

/**
 * @brief Are we writing entropy in the internal energy field ?
 *
 * @return 1 if entropy is in 'internal energy', 0 otherwise.
 */
int writeEntropyFlag() { return 0; }

#endif /* SWIFT_PRESSURE_ENERGY_HYDRO_IO_H */
