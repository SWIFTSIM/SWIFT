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
#ifndef SWIFT_CHEMISTRY_IO_EAGLE_H
#define SWIFT_CHEMISTRY_IO_EAGLE_H

#include "chemistry.h"
#include "io_properties.h"

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
  list[1] =
      io_make_input_field("Metallicity", FLOAT, 1, OPTIONAL, UNIT_CONV_NO_UNITS,
                          parts, chemistry_data.metal_mass_fraction_total);
  list[2] = io_make_input_field("IronMassFracFromSNIa", FLOAT, 1, OPTIONAL,
                                UNIT_CONV_NO_UNITS, parts,
                                chemistry_data.iron_mass_fraction_from_SNIa);

  return 3;
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
  list[0] = io_make_output_field("ElementAbundance", FLOAT,
                                 chemistry_element_count, UNIT_CONV_NO_UNITS,
                                 parts, chemistry_data.metal_mass_fraction);

  list[1] = io_make_output_field(
      "SmoothedElementAbundance", FLOAT, chemistry_element_count,
      UNIT_CONV_NO_UNITS, parts, chemistry_data.smoothed_metal_mass_fraction);

  list[2] =
      io_make_output_field("Metallicity", FLOAT, 1, UNIT_CONV_NO_UNITS, parts,
                           chemistry_data.metal_mass_fraction_total);

  list[3] = io_make_output_field(
      "SmoothedMetallicity", FLOAT, 1, UNIT_CONV_NO_UNITS, parts,
      chemistry_data.smoothed_metal_mass_fraction_total);

  list[4] = io_make_output_field("TotalMassFromSNIa", FLOAT, 1, UNIT_CONV_MASS,
                                 parts, chemistry_data.mass_from_SNIa);

  list[5] = io_make_output_field("MetalMassFracFromSNIa", FLOAT, 1,
                                 UNIT_CONV_NO_UNITS, parts,
                                 chemistry_data.metal_mass_fraction_from_SNIa);

  list[6] = io_make_output_field("TotalMassFromAGB", FLOAT, 1, UNIT_CONV_MASS,
                                 parts, chemistry_data.mass_from_AGB);

  list[7] =
      io_make_output_field("MetalMassFracFromAGB", FLOAT, 1, UNIT_CONV_NO_UNITS,
                           parts, chemistry_data.metal_mass_fraction_from_AGB);

  list[8] = io_make_output_field("TotalMassFromSNII", FLOAT, 1, UNIT_CONV_MASS,
                                 parts, chemistry_data.mass_from_SNII);

  list[9] = io_make_output_field("MetalMassFracFromSNII", FLOAT, 1,
                                 UNIT_CONV_NO_UNITS, parts,
                                 chemistry_data.metal_mass_fraction_from_SNII);

  list[10] =
      io_make_output_field("IronMassFracFromSNIa", FLOAT, 1, UNIT_CONV_NO_UNITS,
                           parts, chemistry_data.iron_mass_fraction_from_SNIa);

  list[11] = io_make_output_field(
      "SmoothedIronMassFracFromSNIa", FLOAT, 1, UNIT_CONV_NO_UNITS, parts,
      chemistry_data.smoothed_iron_mass_fraction_from_SNIa);

  return 12;
}

/**
 * @brief Specifies which star particle fields to write to a dataset
 *
 * @param sparts The star particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
INLINE static int chemistry_write_sparticles(const struct spart* sparts,
                                             struct io_props* list) {

  /* List what we want to write */
  list[0] = io_make_output_field("ElementAbundance", FLOAT,
                                 chemistry_element_count, UNIT_CONV_NO_UNITS,
                                 sparts, chemistry_data.metal_mass_fraction);

  list[1] = io_make_output_field(
      "SmoothedElementAbundance", FLOAT, chemistry_element_count,
      UNIT_CONV_NO_UNITS, sparts, chemistry_data.smoothed_metal_mass_fraction);

  list[2] =
      io_make_output_field("Metallicity", FLOAT, 1, UNIT_CONV_NO_UNITS, sparts,
                           chemistry_data.metal_mass_fraction_total);

  list[3] = io_make_output_field(
      "SmoothedMetallicity", FLOAT, 1, UNIT_CONV_NO_UNITS, sparts,
      chemistry_data.smoothed_metal_mass_fraction_total);

  list[4] = io_make_output_field("TotalMassFromSNIa", FLOAT, 1, UNIT_CONV_MASS,
                                 sparts, chemistry_data.mass_from_SNIa);

  list[5] = io_make_output_field("MetalMassFracFromSNIa", FLOAT, 1,
                                 UNIT_CONV_NO_UNITS, sparts,
                                 chemistry_data.metal_mass_fraction_from_SNIa);

  list[6] = io_make_output_field("TotalMassFromAGB", FLOAT, 1, UNIT_CONV_MASS,
                                 sparts, chemistry_data.mass_from_AGB);

  list[7] =
      io_make_output_field("MetalMassFracFromAGB", FLOAT, 1, UNIT_CONV_NO_UNITS,
                           sparts, chemistry_data.metal_mass_fraction_from_AGB);

  list[8] = io_make_output_field("TotalMassFromSNII", FLOAT, 1, UNIT_CONV_MASS,
                                 sparts, chemistry_data.mass_from_SNII);

  list[9] = io_make_output_field("MetalMassFracFromSNII", FLOAT, 1,
                                 UNIT_CONV_NO_UNITS, sparts,
                                 chemistry_data.metal_mass_fraction_from_SNII);

  list[10] =
      io_make_output_field("IronMassFracFromSNIa", FLOAT, 1, UNIT_CONV_NO_UNITS,
                           sparts, chemistry_data.iron_mass_fraction_from_SNIa);

  list[11] = io_make_output_field(
      "SmoothedIronMassFracFromSNIa", FLOAT, 1, UNIT_CONV_NO_UNITS, sparts,
      chemistry_data.smoothed_iron_mass_fraction_from_SNIa);

  return 12;
}

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grp The HDF5 group in which to write
 */
INLINE static void chemistry_write_flavour(hid_t h_grp) {

  io_write_attribute_s(h_grp, "Chemistry Model", "EAGLE");
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    char buffer[20];
    sprintf(buffer, "Element %d", elem);
    io_write_attribute_s(
        h_grp, buffer,
        chemistry_get_element_name((enum chemistry_element)elem));
  }
}
#endif

#endif /* SWIFT_CHEMISTRY_IO_EAGLE_H */
