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

#include "chemistry_struct.h"
#include "error.h"
#include "io_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

/**
 * @brief Return a string containing the name of a given #chemistry_element.
 */
__attribute__((always_inline)) INLINE static const char*
chemistry_get_element_name(enum chemistry_element elem) {

  static const char* chemistry_element_names[chemistry_element_count] = {
      "Oxygen",    "Magnesium", "Sulfur", "Iron",    "Zinc",
      "Strontium", "Yttrium",   "Barium", "Europium"};

  return chemistry_element_names[elem];
}

/**
 * @brief Specifies which particle fields to read from a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to read.
 *
 * @return Returns the number of fields to read.
 */
INLINE static int chemistry_read_particles(struct part* parts,
                                           struct io_props* list) {

  /* List what we want to read */
  list[0] = io_make_input_field(
      "ElementAbundance", FLOAT, chemistry_element_count, OPTIONAL,
      UNIT_CONV_NO_UNITS, parts, chemistry_data.metal_mass_fraction);
  list[1] = io_make_input_field("Z", FLOAT, 1, OPTIONAL, UNIT_CONV_NO_UNITS,
                                parts, chemistry_data.Z);

  return 2;
}

/**
 * @brief Specifies which particle fields to write to a dataset
 *
 * @param parts The particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
INLINE static int chemistry_write_particles(const struct part* parts,
                                            struct io_props* list) {

  /* List what we want to write */
  list[0] = io_make_output_field(
      "SmoothedElementAbundance", FLOAT, chemistry_element_count,
      UNIT_CONV_NO_UNITS, parts, chemistry_data.smoothed_metal_mass_fraction);

  list[1] = io_make_output_field("Z", FLOAT, 1, UNIT_CONV_NO_UNITS, parts,
                                 chemistry_data.Z);

  list[2] = io_make_output_field("ElementAbundance", FLOAT,
                                 chemistry_element_count, UNIT_CONV_NO_UNITS,
                                 parts, chemistry_data.metal_mass_fraction);

  return 3;
}

/**
 * @brief Specifies which sparticle fields to write to a dataset
 *
 * @param sparts The sparticle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
INLINE static int chemistry_write_sparticles(const struct spart* sparts,
                                             struct io_props* list) {

  /* List what we want to write */
  list[0] = io_make_output_field(
      "SmoothedElementAbundance", FLOAT, chemistry_element_count,
      UNIT_CONV_NO_UNITS, sparts, chemistry_data.smoothed_metal_mass_fraction);

  list[1] = io_make_output_field("Z", FLOAT, 1, UNIT_CONV_NO_UNITS, sparts,
                                 chemistry_data.Z);

  list[2] = io_make_output_field("ElementAbundance", FLOAT,
                                 chemistry_element_count, UNIT_CONV_NO_UNITS,
                                 sparts, chemistry_data.metal_mass_fraction);

  return 3;
}

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grp The HDF5 group in which to write
 */
INLINE static void chemistry_write_flavour(hid_t h_grp) {

  io_write_attribute_s(h_grp, "Chemistry Model", "GEAR");
  for (enum chemistry_element i = chemistry_element_O;
       i < chemistry_element_count; i++) {
    char buffer[20];
    sprintf(buffer, "Element %d", (int)i);
    io_write_attribute_s(h_grp, buffer, chemistry_get_element_name(i));
  }
}
#endif

#endif /* SWIFT_CHEMISTRY_IO_GEAR_H */
