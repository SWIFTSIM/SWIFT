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
#ifndef SWIFT_CHEMISTRY_IO_GEAR_H
#define SWIFT_CHEMISTRY_IO_GEAR_H

#include "chemistry_struct.h"
#include "engine.h"
#include "error.h"
#include "feedback.h"
#include "io_properties.h"
#include "parser.h"
#include "part.h"
#include "physical_constants.h"
#include "units.h"

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
      "MetalMassFraction", DOUBLE, GEAR_CHEMISTRY_ELEMENT_COUNT, OPTIONAL,
      UNIT_CONV_NO_UNITS, parts, chemistry_data.metal_mass);

  return 1;
}

INLINE static void convert_gas_metals(const struct engine* e,
                                      const struct part* p,
                                      const struct xpart* xp, double* ret) {

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    ret[i] = p->chemistry_data.metal_mass[i] / hydro_get_mass(p);
  }
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

  /* List what we want to write */
  list[0] = io_make_output_field(
      "SmoothedMetalMassFractions", DOUBLE, GEAR_CHEMISTRY_ELEMENT_COUNT,
      UNIT_CONV_NO_UNITS, 0.f, parts,
      chemistry_data.smoothed_metal_mass_fraction,
      "Mass fraction of each element smoothed over the neighbors");

  list[1] = io_make_output_field_convert_part(
      "MetalMassFractions", DOUBLE, GEAR_CHEMISTRY_ELEMENT_COUNT,
      UNIT_CONV_NO_UNITS, 0.f, parts, xparts, convert_gas_metals,
      "Mass fraction of each element");

  return 2;
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
      "MetalMassFractions", DOUBLE, GEAR_CHEMISTRY_ELEMENT_COUNT,
      UNIT_CONV_NO_UNITS, 0.f, sparts, chemistry_data.metal_mass_fraction,
      "Mass fraction of each element");

  return 1;
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

  /* No fields to write here */
  return 0;
}

#ifdef HAVE_HDF5

/**
 * @brief Writes the current model of SPH to the file
 * @param h_grp The HDF5 group in which to write
 * @param h_grp_columns The HDF5 group containing named columns
 * @param e The #engine.
 */
INLINE static void chemistry_write_flavour(hid_t h_grp, hid_t h_grp_columns,
                                           const struct engine* e) {
  io_write_attribute_s(h_grp, "Chemistry model", "GEAR");
  io_write_attribute_d(h_grp, "Chemistry element count",
                       GEAR_CHEMISTRY_ELEMENT_COUNT);
#if defined(FEEDBACK_GEAR) || FEEDBACK_GEAR_MECHANICAL_MODE >= 1
  /* Without feedback, the elements are meaningless */
  const int with_feedback = e->policy & engine_policy_feedback;
  if (!with_feedback) return;

  const char* element_names = e->feedback_props->stellar_model.elements_name;

  /* Add to the named columns */
  hsize_t dims[1] = {GEAR_CHEMISTRY_ELEMENT_COUNT};
  hid_t type = H5Tcopy(H5T_C_S1);
  H5Tset_size(type, GEAR_LABELS_SIZE);
  hid_t space = H5Screate_simple(1, dims, NULL);
  hid_t dset = H5Dcreate(h_grp_columns, "MetalMassFractions", type, space,
                         H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, element_names);
  H5Dclose(dset);
  dset = H5Dcreate(h_grp_columns, "SmoothedMetalMassFractions", type, space,
                   H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dset, type, H5S_ALL, H5S_ALL, H5P_DEFAULT, element_names);
  H5Dclose(dset);

  H5Tclose(type);

  /* Write the solar abundances and the elements */
  /* Create the group */
  hid_t h_sol_ab = H5Gcreate(h_grp, "SolarAbundances", H5P_DEFAULT, H5P_DEFAULT,
                             H5P_DEFAULT);
  if (h_sol_ab < 0) error("Error while creating the SolarAbundances group\n");

  /* Write all the elements as attributes */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    const char* name = stellar_evolution_get_element_name(
        &e->feedback_props->stellar_model, i);

    io_write_attribute_f(h_sol_ab, name, e->chemistry->solar_abundances[i]);
  }

  /* Close group */
  H5Gclose(h_sol_ab);

#endif
}
#endif

#endif /* SWIFT_CHEMISTRY_IO_GEAR_H */
