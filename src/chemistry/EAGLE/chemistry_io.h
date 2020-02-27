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
  list[0] = io_make_output_field(
      "ElementMassFractions", FLOAT, chemistry_element_count,
      UNIT_CONV_NO_UNITS, 0.f, parts, chemistry_data.metal_mass_fraction,
      "Fractions of the particles' masses that are in the given element");

  list[1] = io_make_output_field(
      "SmoothedElementMassFractions", FLOAT, chemistry_element_count,
      UNIT_CONV_NO_UNITS, 0.f, parts,
      chemistry_data.smoothed_metal_mass_fraction,
      "Smoothed fractions of the particles' masses that are "
      "in the given element");

  list[2] = io_make_output_field(
      "MetalMassFractions", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, parts,
      chemistry_data.metal_mass_fraction_total,
      "Fractions of the particles' masses that are in metals");

  list[3] = io_make_output_field(
      "SmoothedMetalMassFractions", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, parts,
      chemistry_data.smoothed_metal_mass_fraction_total,
      "Smoothed fractions of the particles masses that are in metals");

  list[4] = io_make_output_field(
      "MassesFromSNIa", FLOAT, 1, UNIT_CONV_MASS, 0.f, parts,
      chemistry_data.mass_from_SNIa,
      "Masses of gas that have been produced by SNIa stars");

  list[5] = io_make_output_field("MetalMassFractionsFromSNIa", FLOAT, 1,
                                 UNIT_CONV_NO_UNITS, 0.f, parts,
                                 chemistry_data.metal_mass_fraction_from_SNIa,
                                 "Fractions of the particles' masses that are "
                                 "in metals produced by SNIa stars");

  list[6] = io_make_output_field(
      "MassesFromAGB", FLOAT, 1, UNIT_CONV_MASS, 0.f, parts,
      chemistry_data.mass_from_AGB,
      "Masses of gas that have been produced by AGN stars");

  list[7] = io_make_output_field("MetalMassFractionsFromAGB", FLOAT, 1,
                                 UNIT_CONV_NO_UNITS, 0., parts,
                                 chemistry_data.metal_mass_fraction_from_AGB,
                                 "Fractions of the particles' masses that are "
                                 "in metals produced by AGB stars");

  list[8] = io_make_output_field(
      "MassesFromSNII", FLOAT, 1, UNIT_CONV_MASS, 0.f, parts,
      chemistry_data.mass_from_SNII,
      "Masses of gas that have been produced by SNII stars");

  list[9] = io_make_output_field("MetalMassFractionsFromSNII", FLOAT, 1,
                                 UNIT_CONV_NO_UNITS, 0.f, parts,
                                 chemistry_data.metal_mass_fraction_from_SNII,
                                 "Fractions of the particles' masses that are "
                                 "in metals produced by SNII stars");

  list[10] = io_make_output_field("IronMassFractionsFromSNIa", FLOAT, 1,
                                  UNIT_CONV_NO_UNITS, 0.f, parts,
                                  chemistry_data.iron_mass_fraction_from_SNIa,
                                  "Fractions of the particles' masses that are "
                                  "in iron produced by SNIa stars");

  list[11] = io_make_output_field(
      "SmoothedIronMassFractionsFromSNIa", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
      parts, chemistry_data.smoothed_iron_mass_fraction_from_SNIa,
      "Smoothed fractions of the particles' masses that are "
      "in iron produced by SNIa stars");

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
  list[0] = io_make_output_field(
      "ElementMassFractions", FLOAT, chemistry_element_count,
      UNIT_CONV_NO_UNITS, 0.f, sparts, chemistry_data.metal_mass_fraction,
      "Fractions of the particles' masses that are in the given element");

  list[1] = io_make_output_field(
      "SmoothedElementMassFractions", FLOAT, chemistry_element_count,
      UNIT_CONV_NO_UNITS, 0.f, sparts,
      chemistry_data.smoothed_metal_mass_fraction,
      "Smoothed fractions of the particles' masses that are "
      "in the given element");

  list[2] = io_make_output_field(
      "MetalMassFractions", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      chemistry_data.metal_mass_fraction_total,
      "Fractions of the particles' masses that are in metals");

  list[3] = io_make_output_field(
      "SmoothedMetalMassFractions", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      chemistry_data.smoothed_metal_mass_fraction_total,
      "Smoothed fractions of the particles masses that are in metals");

  list[4] = io_make_output_field(
      "MassesFromSNIa", FLOAT, 1, UNIT_CONV_MASS, 0.f, sparts,
      chemistry_data.mass_from_SNIa,
      "Masses of gas that have been produced by SNIa stars");

  list[5] = io_make_output_field("MetalMassFractionsFromSNIa", FLOAT, 1,
                                 UNIT_CONV_NO_UNITS, 0.f, sparts,
                                 chemistry_data.metal_mass_fraction_from_SNIa,
                                 "Fractions of the particles' masses that are "
                                 "in metals produced by SNIa stars");

  list[6] = io_make_output_field(
      "MassesFromAGB", FLOAT, 1, UNIT_CONV_MASS, 0.f, sparts,
      chemistry_data.mass_from_AGB,
      "Masses of gas that have been produced by AGN stars");

  list[7] = io_make_output_field("MetalMassFractionsFromAGB", FLOAT, 1,
                                 UNIT_CONV_NO_UNITS, 0., sparts,
                                 chemistry_data.metal_mass_fraction_from_AGB,
                                 "Fractions of the particles' masses that are "
                                 "in metals produced by AGB stars");

  list[8] = io_make_output_field(
      "MassesFromSNII", FLOAT, 1, UNIT_CONV_MASS, 0.f, sparts,
      chemistry_data.mass_from_SNII,
      "Masses of gas that have been produced by SNII stars");

  list[9] = io_make_output_field("MetalMassFractionsFromSNII", FLOAT, 1,
                                 UNIT_CONV_NO_UNITS, 0.f, sparts,
                                 chemistry_data.metal_mass_fraction_from_SNII,
                                 "Fractions of the particles' masses that are "
                                 "in metals produced by SNII stars");

  list[10] = io_make_output_field("IronMassFractionsFromSNIa", FLOAT, 1,
                                  UNIT_CONV_NO_UNITS, 0.f, sparts,
                                  chemistry_data.iron_mass_fraction_from_SNIa,
                                  "Fractions of the particles' masses that are "
                                  "in iron produced by SNIa stars");

  list[11] = io_make_output_field(
      "SmoothedIronMassFractionsFromSNIa", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f,
      sparts, chemistry_data.smoothed_iron_mass_fraction_from_SNIa,
      "Smoothed fractions of the particles' masses that are "
      "in iron produced by SNIa stars");

  return 12;
}

/**
 * @brief Specifies which black hole particle fields to write to a dataset
 *
 * @param bparts The black hole particle array.
 * @param list The list of i/o properties to write.
 *
 * @return Returns the number of fields to write.
 */
INLINE static int chemistry_write_bparticles(const struct bpart* bparts,
                                             struct io_props* list) {

  /* List what we want to write */
  list[0] = io_make_output_field(
      "ElementMasses", FLOAT, chemistry_element_count, UNIT_CONV_MASS, 0.f,
      bparts, chemistry_data.metal_mass,
      "Masses of the BH particles in a given element");

  list[1] = io_make_output_field("MetalMasses", FLOAT, chemistry_element_count,
                                 UNIT_CONV_MASS, 0.f, bparts,
                                 chemistry_data.metal_mass_total,
                                 "Masses of the BH particles in a metals");

  list[2] = io_make_output_field(
      "MassesFromSNIa", FLOAT, 1, UNIT_CONV_MASS, 0.f, bparts,
      chemistry_data.mass_from_SNIa,
      "Masses of the BH particles that have been produced by SNIa stars");

  list[3] = io_make_output_field(
      "MassesFromSNII", FLOAT, 1, UNIT_CONV_MASS, 0.f, bparts,
      chemistry_data.mass_from_SNII,
      "Masses of the BH particles that have been produced by SNII stars");

  list[4] = io_make_output_field(
      "MassesFromAGB", FLOAT, 1, UNIT_CONV_MASS, 0.f, bparts,
      chemistry_data.mass_from_AGB,
      "Masses of the BH particles that have been produced by AGB stars");

  list[5] =
      io_make_output_field("MetalMassesFromSNIa", FLOAT, 1, UNIT_CONV_MASS, 0.f,
                           bparts, chemistry_data.metal_mass_from_SNIa,
                           "Masses of the BH particles in metals that have "
                           "been produced by SNIa stars");

  list[6] =
      io_make_output_field("MetalMassesFromSNII", FLOAT, 1, UNIT_CONV_MASS, 0.f,
                           bparts, chemistry_data.metal_mass_from_SNII,
                           "Masses of the BH particles in metals that have "
                           "been produced by SNII stars");

  list[7] =
      io_make_output_field("MetalMassesFromAGB", FLOAT, 1, UNIT_CONV_MASS, 0.f,
                           bparts, chemistry_data.metal_mass_from_AGB,
                           "Masses of the BH particles in metals that have "
                           "been produced by AGB stars");

  list[8] =
      io_make_output_field("IronMassesFromSNIa", FLOAT, 1, UNIT_CONV_MASS, 0.f,
                           bparts, chemistry_data.iron_mass_from_SNIa,
                           "Masses of the BH particles in iron that have been "
                           "produced by SNIa stars");

  return 9;
}

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grp The HDF5 group in which to write
 */
INLINE static void chemistry_write_flavour(hid_t h_grp) {

  io_write_attribute_s(h_grp, "Chemistry Model", "EAGLE");
  io_write_attribute_d(h_grp, "Chemistry element count",
                       chemistry_element_count);
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
