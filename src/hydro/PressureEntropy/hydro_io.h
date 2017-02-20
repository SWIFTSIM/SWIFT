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
#ifndef SWIFT_PRESSURE_ENTROPY_HYDRO_IO_H
#define SWIFT_PRESSURE_ENTROPY_HYDRO_IO_H

/**
 * @file PressureEntropy/hydro_io.h
 * @brief Pressure-Entropy implementation of SPH (i/o routines)
 *
 * The thermal variable is the entropy (S) and the entropy is smoothed over
 * contact discontinuities to prevent spurious surface tension.
 *
 * Follows eqautions (19), (21) and (22) of Hopkins, P., MNRAS, 2013,
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
  list[4] =
      io_make_input_field("InternalEnergy", FLOAT, 1, COMPULSORY,
                          UNIT_CONV_ENTROPY_PER_UNIT_MASS, parts, entropy);
  list[5] = io_make_input_field("ParticleIDs", ULONGLONG, 1, COMPULSORY,
                                UNIT_CONV_NO_UNITS, parts, id);
  list[6] = io_make_input_field("Accelerations", FLOAT, 3, OPTIONAL,
                                UNIT_CONV_ACCELERATION, parts, a_hydro);
  list[7] = io_make_input_field("Density", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_DENSITY, parts, rho);
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

  *num_fields = 11;

  /* List what we want to write */
  list[0] = io_make_output_field("Coordinates", DOUBLE, 3, UNIT_CONV_LENGTH,
                                 parts, x);
  list[1] =
      io_make_output_field("Velocities", FLOAT, 3, UNIT_CONV_SPEED, parts, v);
  list[2] =
      io_make_output_field("Masses", FLOAT, 1, UNIT_CONV_MASS, parts, mass);
  list[3] = io_make_output_field("SmoothingLength", FLOAT, 1, UNIT_CONV_LENGTH,
                                 parts, h);
  list[4] = io_make_output_field_convert_part("InternalEnergy", FLOAT, 1,
                                              UNIT_CONV_ENERGY_PER_UNIT_MASS,
                                              parts, entropy, convert_u);
  list[5] = io_make_output_field("ParticleIDs", ULONGLONG, 1,
                                 UNIT_CONV_NO_UNITS, parts, id);
  list[6] = io_make_output_field("Acceleration", FLOAT, 3,
                                 UNIT_CONV_ACCELERATION, parts, a_hydro);
  list[7] =
      io_make_output_field("Density", FLOAT, 1, UNIT_CONV_DENSITY, parts, rho);
  list[8] = io_make_output_field(
      "Entropy", FLOAT, 1, UNIT_CONV_ENTROPY_PER_UNIT_MASS, parts, entropy);
  list[9] = io_make_output_field_convert_part(
      "Pressure", FLOAT, 1, UNIT_CONV_PRESSURE, parts, rho, convert_P);
  list[10] = io_make_output_field("WeightedDensity", FLOAT, 1,
                                  UNIT_CONV_DENSITY, parts, rho_bar);
}

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grpsph The HDF5 group in which to write
 */
void writeSPHflavour(hid_t h_grpsph) {

  /* Viscosity and thermal conduction */
  /* Nothing in this minimal model... */
  io_write_attribute_s(h_grpsph, "Thermal Conductivity Model", "No treatment");
  io_write_attribute_s(
      h_grpsph, "Viscosity Model",
      "as in Springel (2005), i.e. Monaghan (1992) with Balsara (1995) switch");
  io_write_attribute_f(h_grpsph, "Viscosity alpha", const_viscosity_alpha);
  io_write_attribute_f(h_grpsph, "Viscosity beta", 3.f);

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

#endif /* SWIFT_PRESSURE_ENTROPY_HYDRO_IO_H */
