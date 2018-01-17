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
#ifndef SWIFT_CHEMISTRY_IO_GEAR_H
#define SWIFT_CHEMISTRY_IO_GEAR_H

#include "io_properties.h"

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 *
 * @return Returns the number of fields to read.
 */
int chemistry_read_particles(struct part* parts, struct io_props* list) {

  /* List what we want to read */
  list[0] =
      io_make_input_field("HeDensity", FLOAT, 1, COMPULSORY, UNIT_CONV_DENSITY,
                          parts, chemistry_data.he_density);

  return 1;
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
int chemistry_write_particles(const struct part* parts, struct io_props* list) {

  /* List what we want to write */
  list[0] = io_make_output_field("HeDensity", FLOAT, 1, UNIT_CONV_DENSITY,
                                 parts, chemistry_data.he_density);

  return 1;
}

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grpsph The HDF5 group in which to write
 */
void chemistry_write_flavour(hid_t h_grpsph) {

  io_write_attribute_s(h_grpsph, "Chemistry Model", "GEAR");
}

#endif /* SWIFT_CHEMISTRY_IO_GEAR_H */
