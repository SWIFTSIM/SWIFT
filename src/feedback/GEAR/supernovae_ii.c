/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2019 Loic Hausammann (loic.hausammann@epfl.ch)
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

/* Include header */
#include "supernovae_ii.h"

/* Local headers */
#include "hdf5_functions.h"
#include "interpolation.h"
#include "stellar_evolution.h"
#include "stellar_evolution_struct.h"

/**
 * @brief Print the supernovae II model.
 *
 * @param snii The #supernovae_ii.
 */
void supernovae_ii_print(const struct supernovae_ii *snii) {

  /* Only the master print */
  if (engine_rank != 0) {
    return;
  }

  message("Mass range for SNII = [%g, %g]", snii->mass_min, snii->mass_max);
}

/**
 * @brief Check if the given mass is able to produce a SNII.
 *
 * @param snii The #supernovae_ii model.
 * @param m_low The lower mass.
 * @param m_high The higher mass
 *
 * @return If the mass is in the range of SNII.
 */
int supernovae_ii_can_explode(const struct supernovae_ii *snii, float m_low,
                              float m_high) {

  if (m_high < snii->mass_min || m_low > snii->mass_max) return 0;

  return 1;
}

/**
 * @brief Compute the number of supernovae II per unit of mass (equation 3.47 in
 * Poirier 2004).
 *
 * @param snii The #supernovae_ii model.
 * @param m1 The lower mass limit.
 * @param m2 The upper mass limit.
 *
 * @return The number of supernovae II per unit of mass.
 */
float supernovae_ii_get_number_per_unit_mass(const struct supernovae_ii *snii, float m1,
					     float m2) {
#ifdef SWIFT_DEBUG_CHECKS
  if (m1 > m2) error("Mass 1 larger than mass 2 %g > %g.", m1, m2);
#endif

  /* Can we explode SNII? */
  if (!supernovae_ii_can_explode(snii, m1, m2)) {
    return 0.;
  }

  const float mass_min = max(m1, snii->mass_min);
  const float mass_max = min(m2, snii->mass_max);

  const float pow_mass =
      pow(mass_max, snii->exponent) - pow(mass_min, snii->exponent);

  return snii->coef_exp * pow_mass;
};

/**
 * @brief Get the SNII yields per mass (Poirier version).
 *
 * @param snii The #supernovae_ii model.
 * @param m1 The lower mass in log.
 * @param m2 The upper mass in log.
 * @param yields The elements ejected (needs to be allocated).
 */
void supernovae_ii_get_yields(const struct supernovae_ii *snii, float log_m1,
                              float log_m2, float *yields) {

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    float yields_1 = interpolate_1d(&snii->integrated_yields[i], log_m1);
    float yields_2 = interpolate_1d(&snii->integrated_yields[i], log_m2);

    yields[i] = yields_2 - yields_1;
  }
};

/**
 * @brief Get the ejected mass per mass unit.
 *
 * @param snii The #supernovae_ii model.
 * @param m1 The lower mass in log.
 * @param m2 The upper mass in log.
 *
 * @return mass_ejected_processed The mass of processsed elements.
 */
float supernovae_ii_get_ejected_mass_fraction(const struct supernovae_ii *snii,
                                              float log_m1, float log_m2) {

  float mass_ejected_1 = interpolate_1d(&snii->integrated_ejected_mass, log_m1);
  float mass_ejected_2 = interpolate_1d(&snii->integrated_ejected_mass, log_m2);

  return mass_ejected_2 - mass_ejected_1;
};

/**
 * @brief Get the ejected mass (processed) per mass.
 *
 * @param snii The #supernovae_ii model.
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 *
 * @return mass_ejected The mass of non processsed elements.
 */
float supernovae_ii_get_ejected_mass_fraction_processed(
    const struct supernovae_ii *snii, float log_m1, float log_m2) {

  float mass_ejected_1 =
      interpolate_1d(&snii->integrated_ejected_mass_processed, log_m1);
  float mass_ejected_2 =
      interpolate_1d(&snii->integrated_ejected_mass_processed, log_m2);

  return mass_ejected_2 - mass_ejected_1;
};

/**
 * @brief Read an array of SNII yields from the table.
 *
 * @param snii The #supernovae_ii model.
 * @param sm The #stellar_model.
 */
void supernovae_ii_read_yields_array(
    struct supernovae_ii *snii, struct interpolation_1d *interp,
    const struct phys_const *phys_const, const struct stellar_model *sm,
    hid_t group_id, const char *hdf5_dataset_name, hsize_t *previous_count,
    int interpolation_size) {

  /* Now let's get the number of elements */
  /* Open attribute */
  const hid_t h_dataset = H5Dopen(group_id, hdf5_dataset_name, H5P_DEFAULT);
  if (h_dataset < 0)
    error("Error while opening attribute '%s'", hdf5_dataset_name);

  /* Get the number of elements */
  hsize_t count = io_get_number_element_in_dataset(h_dataset);

  /* Check that all the arrays have the same size */
  if (*previous_count != 0 && count != *previous_count) {
    error("The code is not able to deal with yields arrays of different size");
  }
  *previous_count = count;

  /* Read the minimal mass (in log) */
  float log_mass_min = 0;
  io_read_attribute(h_dataset, "min", FLOAT, &log_mass_min);

  /* Read the step size (log step) */
  float step_size = 0;
  io_read_attribute(h_dataset, "step", FLOAT, &step_size);

  /* Close the attribute */
  H5Dclose(h_dataset);

  /* Allocate the memory */
  float *data = (float *)malloc(sizeof(float) * count);
  if (data == NULL)
    error("Failed to allocate the SNII yields for %s.", hdf5_dataset_name);

  /* Read the dataset */
  io_read_array_dataset(group_id, hdf5_dataset_name, FLOAT, data, count);

  /* Integrate the yields */
  initial_mass_function_integrate(&sm->imf, data, count, log_mass_min, step_size);
  // TODO: decrease count in order to keep the same distance between points

  /* Initialize the interpolation */
  interpolate_1d_init(interp, log10(snii->mass_min), log10(snii->mass_max),
                      interpolation_size, log_mass_min, step_size, count, data,
		      boundary_condition_zero_const);


  /* Cleanup the memory */
  free(data);
}

/**
 * @brief Read the SNII yields from the table.
 *
 * The tables are in internal units at the end of this function.
 *
 * @param snii The #supernovae_ii model.
 * @param params The simulation parameters.
 * @param sm The #stellar_model.
 */
void supernovae_ii_read_yields(struct supernovae_ii *snii,
                               struct swift_params *params,
                               const struct phys_const *phys_const,
                               const struct stellar_model *sm) {

  hid_t file_id, group_id;

  hsize_t previous_count = 0;

  const int interpolation_size = parser_get_opt_param_int(
      params, "GEARSupernovaeII:interpolation_size", 200);

  /* Open IMF group */
  h5_open_group(params, "Data/SNII", &file_id, &group_id);

  /* Do all the elements */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {

    /* Get the element name */
    const char *name = stellar_evolution_get_element_name(sm, i);

    /* Read the array */
    supernovae_ii_read_yields_array(snii, &snii->integrated_yields[i],
                                    phys_const, sm, group_id, name,
                                    &previous_count, interpolation_size);
  }

  /* Read the mass ejected */
  supernovae_ii_read_yields_array(
      snii, &snii->integrated_ejected_mass_processed, phys_const, sm, group_id,
      "Ej", &previous_count, interpolation_size);

  /* Read the mass ejected of non processed gas */
  supernovae_ii_read_yields_array(snii, &snii->integrated_ejected_mass,
                                  phys_const, sm, group_id, "Ejnp",
                                  &previous_count, interpolation_size);

  /* Cleanup everything */
  h5_close_group(file_id, group_id);
};

/**
 * @brief Reads the supernovae II parameters from parameters file.
 *
 * @param snii The #supernovae_ii model.
 * @param params The simulation parameters.
 */
void supernovae_ii_read_from_params(struct supernovae_ii *snii,
                                    struct swift_params *params) {

  /* Read the minimal mass of a supernovae */
  snii->mass_min = parser_get_opt_param_float(
      params, "GEARSupernovaeII:min_mass", snii->mass_min);

  /* Read the maximal mass of a supernovae */
  snii->mass_max = parser_get_opt_param_float(
      params, "GEARSupernovaeII:max_mass", snii->mass_max);
}

/**
 * @brief Reads the supernovae II parameters from tables.
 *
 * @param snii The #supernovae_ii model.
 * @param params The simulation parameters.
 */
void supernovae_ii_read_from_tables(struct supernovae_ii *snii,
                                    struct swift_params *params) {

  hid_t file_id, group_id;

  /* Open IMF group */
  h5_open_group(params, "Data/SNII", &file_id, &group_id);

  /* Read the minimal mass of a supernovae */
  io_read_attribute(group_id, "Mmin", FLOAT, &snii->mass_min);

  /* Read the maximal mass of a supernovae */
  io_read_attribute(group_id, "Mmax", FLOAT, &snii->mass_max);

  /* Cleanup everything */
  h5_close_group(file_id, group_id);
}

/**
 * @brief Initialize the #supernovae_ii structure.
 *
 * @param snii The #supernovae_ii model.
 * @param phys_const The #phys_const.
 * @param us The #unit_system.
 * @param params The simulation parameters.
 * @param sm The #stellar_model.
 */
void supernovae_ii_init(struct supernovae_ii *snii,
                        const struct phys_const *phys_const,
                        const struct unit_system *us,
                        struct swift_params *params,
                        const struct stellar_model *sm) {

  /* Read the parameters from the tables */
  supernovae_ii_read_from_tables(snii, params);

  /* Read the parameters from the params file */
  supernovae_ii_read_from_tables(snii, params);

  /* Read the supernovae yields (and apply the units) */
  supernovae_ii_read_yields(snii, params, phys_const, sm);

  /* Get the IMF parameters */
  snii->exponent = initial_mass_function_get_exponent(&sm->imf, snii->mass_min,
                                                      snii->mass_max);
  snii->coef_exp = initial_mass_function_get_coefficient(
      &sm->imf, snii->mass_min, snii->mass_max);
  snii->coef_exp /= snii->exponent;
}

/**
 * @brief Write a supernovae_ii struct to the given FILE as a stream of bytes.
 *
 * Here we are only writing the arrays, everything else has been copied in the
 * feedback.
 *
 * @param snii the struct
 * @param stream the file stream
 */
void supernovae_ii_dump(const struct supernovae_ii *snii, FILE *stream,
                        const struct stellar_model *sm) {

  /* Dump the yields. */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    if (snii->integrated_yields[i].data == NULL) {
      continue;
    }

    const char *name = stellar_evolution_get_element_name(sm, i);
    restart_write_blocks((void *)snii->integrated_yields[i].data, sizeof(float),
                         snii->integrated_yields[i].N, stream, name, name);
  }

  /*! Dump the processed mass. */
  if (snii->integrated_ejected_mass_processed.data != NULL) {
    restart_write_blocks((void *)snii->integrated_ejected_mass_processed.data,
                         sizeof(float),
                         snii->integrated_ejected_mass_processed.N, stream,
                         "processed_mass", "processed_mass");
  }

  /*! Dump the non processed mass. */
  if (snii->integrated_ejected_mass.data != NULL) {
    restart_write_blocks((void *)snii->integrated_ejected_mass.data,
                         sizeof(float), snii->integrated_ejected_mass.N, stream,
                         "non_processed_mass", "non_processed_mass");
  }
}

/**
 * @brief Restore a supernovae_ii struct from the given FILE as a stream of
 * bytes.
 *
 * Here we are only writing the arrays, everything else has been copied in the
 * feedback.
 *
 * @param snii the struct
 * @param stream the file stream
 */
void supernovae_ii_restore(struct supernovae_ii *snii, FILE *stream,
                           const struct stellar_model *sm) {

  /* Restore the yields */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    if (snii->integrated_yields[i].data == NULL) {
      continue;
    }

    const char *name = stellar_evolution_get_element_name(sm, i);
    snii->integrated_yields[i].data =
        (float *)malloc(sizeof(float) * snii->integrated_yields[i].N);

    restart_read_blocks((void *)snii->integrated_yields[i].data, sizeof(float),
                        snii->integrated_yields[i].N, stream, NULL, name);
  }

  /* Restore the processed mass */

  if (snii->integrated_ejected_mass_processed.data != NULL) {
    snii->integrated_ejected_mass_processed.data = (float *)malloc(
        sizeof(float) * snii->integrated_ejected_mass_processed.N);

    restart_read_blocks((void *)snii->integrated_ejected_mass_processed.data,
                        sizeof(float),
                        snii->integrated_ejected_mass_processed.N, stream, NULL,
                        "processed_mass");
  }

  /* Restore the non processed mass */
  if (snii->integrated_ejected_mass.data != NULL) {
    snii->integrated_ejected_mass.data =
        (float *)malloc(sizeof(float) * snii->integrated_ejected_mass.N);

    restart_read_blocks((void *)snii->integrated_ejected_mass.data,
                        sizeof(float), snii->integrated_ejected_mass.N, stream,
                        NULL, "non_processed_mass");
  }
}

/**
 * @brief Clean the allocated memory.
 *
 * @param snii the #supernovae_ii.
 */
void supernovae_ii_clean(struct supernovae_ii *snii) {

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    interpolate_1d_free(&snii->integrated_yields[i]);
  }

  interpolate_1d_free(&snii->integrated_ejected_mass_processed);
  interpolate_1d_free(&snii->integrated_ejected_mass);
}
