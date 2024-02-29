/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_CHEMISTRY_IO_SIMBA_H
#define SWIFT_CHEMISTRY_IO_SIMBA_H

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
 * @param xparts The extra particle array.
 * @param list The list of i/o properties to write.
 * @param with_cosmology Are we running with cosmology?
 *
 * @return Returns the number of fields to write.
 */
INLINE static int chemistry_write_particles(const struct part* parts,
                                            const struct xpart* xparts,
                                            struct io_props* list,
                                            const int with_cosmology) {

  int num=0;

  /* List what we want to write */
  list[num] = io_make_output_field(
      "ElementMassFractions", FLOAT, chemistry_element_count,
      UNIT_CONV_NO_UNITS, 0.f, parts, chemistry_data.metal_mass_fraction,
      "Fractions of the particles' masses that are in the given element");
  num++;

  list[num] = io_make_output_field(
      "MetalMassFractions", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, parts,
      chemistry_data.metal_mass_fraction_total,
      "Fractions of the particles' masses that are in metals");
  num++;

  return num;
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

  int num=0;

  /* List what we want to write */
  list[num] = io_make_output_field(
      "ElementMassFractions", FLOAT, chemistry_element_count,
      UNIT_CONV_NO_UNITS, 0.f, sparts, chemistry_data.metal_mass_fraction,
      "Fractions of the particles' masses that are in the given element");

  list[num] = io_make_output_field(
      "MetalMassFractions", FLOAT, 1, UNIT_CONV_NO_UNITS, 0.f, sparts,
      chemistry_data.metal_mass_fraction_total,
      "Fractions of the particles' masses that are in metals");

  return num;
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

  int num=0;

  /* List what we want to write */
  list[num] = io_make_output_field(
      "ElementMasses", FLOAT, chemistry_element_count, UNIT_CONV_MASS, 0.f,
      bparts, chemistry_data.metal_mass,
      "Masses of the BH particles in a given element");
  num++;

  list[num] = io_make_output_field("MetalMasses", FLOAT, chemistry_element_count,
                                 UNIT_CONV_MASS, 0.f, bparts,
                                 chemistry_data.metal_mass_total,
                                 "Masses of the BH particles in a metals");
  num++;

  return num;
}

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of chemistry to the file
 * @param h_grp The HDF5 group in which to write
 * @param h_grp_columns The HDF5 group containing named columns
 * @param e The #engine.
 */
INLINE static void chemistry_write_flavour(hid_t h_grp, hid_t h_grp_columns,
                                           const struct engine* e) {

  /* Write the chemistry model */
  io_write_attribute_s(h_grp, "Chemistry Model", "SIMBA");

  /* Create an array of element names */
  const int element_name_length = 32;
  char element_names[chemistry_element_count][element_name_length];
  for (int elem = 0; elem < chemistry_element_count; ++elem) {
    sprintf(element_names[elem], "%s",
            chemistry_get_element_name((enum chemistry_element)elem));
  }

  /* Add to the named columns */
  hsize_t dims[1] = {chemistry_element_count};
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, element_name_length);
  hid_t space = H5Screate_simple(1, dims, NULL);
  hid_t dset = H5Dcreate(h_grp_columns, "ElementMassFractions", type, space,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, element_names[0]);
  H5Dclose(dset);
  dset = H5Dcreate(h_grp_columns, "SmoothedElementMassFractions", type, space,
                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, element_names[0]);
  H5Dclose(dset);

  H5Tclose(type);
  H5Sclose(space);
}
#endif

#endif /* SWIFT_CHEMISTRY_IO_SIMBA_H */
