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
#ifndef SWIFT_DEFAULT_HYDRO_IO_H
#define SWIFT_DEFAULT_HYDRO_IO_H

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

  *num_fields = 8;

  /* List what we want to read */
  list[0] = io_make_input_field("Coordinates", DOUBLE, 3, COMPULSORY,
                                UNIT_CONV_LENGTH, parts, x);
  list[1] = io_make_input_field("Velocities", FLOAT, 3, COMPULSORY,
                                UNIT_CONV_SPEED, parts, v);
  list[2] = io_make_input_field("Masses", FLOAT, 1, COMPULSORY, UNIT_CONV_MASS,
                                parts, mass);
  list[3] = io_make_input_field("SmoothingLength", FLOAT, 1, COMPULSORY,
                                UNIT_CONV_LENGTH, parts, h);
  list[4] = io_make_input_field("InternalEnergy", FLOAT, 1, COMPULSORY,
                                UNIT_CONV_ENERGY_PER_UNIT_MASS, parts, u);
  list[5] = io_make_input_field("ParticleIDs", ULONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, parts, id);
  list[6] = io_make_input_field("Accelerations", FLOAT, 3, OPTIONAL,
                                UNIT_CONV_ACCELERATION, parts, a_hydro);
  list[7] = io_make_input_field("Density", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_DENSITY, parts, rho);
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
  list[2] =
      io_make_output_field("Masses", FLOAT, 1, UNIT_CONV_MASS, parts, mass);
  list[3] = io_make_output_field("SmoothingLength", FLOAT, 1, UNIT_CONV_LENGTH,
                                 parts, h);
  list[4] = io_make_output_field("InternalEnergy", FLOAT, 1,
                                 UNIT_CONV_ENERGY_PER_UNIT_MASS, parts, u);
  list[5] = io_make_output_field("ParticleIDs", ULONGLONG, 1,
                                 UNIT_CONV_NO_UNITS, parts, id);
  list[6] = io_make_output_field("Acceleration", FLOAT, 3,
                                 UNIT_CONV_ACCELERATION, parts, a_hydro);
  list[7] =
      io_make_output_field("Density", FLOAT, 1, UNIT_CONV_DENSITY, parts, rho);
}

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grpsph The HDF5 group in which to write
 */
void writeSPHflavour(hid_t h_grpsph) {

  /* Viscosity and thermal conduction */
  io_write_attribute_s(h_grpsph, "Thermal Conductivity Model",
                       "Price (2008) without switch");
  io_write_attribute_f(h_grpsph, "Thermal Conductivity alpha",
                       const_conductivity_alpha);
  io_write_attribute_s(
      h_grpsph, "Viscosity Model",
      "Morris & Monaghan (1997), Rosswog, Davies, Thielemann & "
      "Piran (2000) with additional Balsara (1995) switch");
  io_write_attribute_f(h_grpsph, "Viscosity alpha_min",
                       const_viscosity_alpha_min);
  io_write_attribute_f(h_grpsph, "Viscosity alpha_max",
                       const_viscosity_alpha_max);
  io_write_attribute_f(h_grpsph, "Viscosity beta", 2.f);
  io_write_attribute_f(h_grpsph, "Viscosity decay length",
                       const_viscosity_length);

  /* Time integration properties */
  io_write_attribute_f(h_grpsph, "Maximal Delta u change over dt",
                       const_max_u_change);
}

/**
 * @brief Are we writing entropy in the internal energy field ?
 *
 * @return 1 if entropy is in 'internal energy', 0 otherwise.
 */
int writeEntropyFlag() { return 0; }

#endif /* SWIFT_DEFAULT_HYDRO_IO_H */
