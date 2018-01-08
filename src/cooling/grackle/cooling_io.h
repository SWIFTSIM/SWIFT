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

#include "io_properties.h"

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 * @param num_fields The number of i/o fields to read.
 */
void cooling_read_particles(struct part* parts, struct io_props* list,
                          int* num_fields) {

  list += *num_fields;
  *num_fields += 1;

  /* List what we want to read */
  list[0] = io_make_input_field("HeDensity", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_DENSITY, parts, cooling_data.He_density);

}


/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to write.
 * @param num_fields The number of i/o fields to write.
 */
void cooling_write_particles(const struct part* parts, struct io_props* list,
                           int* num_fields) {

  int i = *num_fields;
  /* List what we want to write */
  list[i] = io_make_output_field(
      "HeDensity", FLOAT, 1, UNIT_CONV_DENSITY, parts, cooling_data.He_density);
  i++;

  *num_fields = i;
}


/**
 * @brief Writes the current model of SPH to the file
 * @param h_grpsph The HDF5 group in which to write
 */
void writeCoolingFlavor(hid_t h_grpsph) {

  /* Viscosity and thermal conduction */
  io_write_attribute_s(
      h_grpsph, "Chemistry Model",
      "Grackle");
}
