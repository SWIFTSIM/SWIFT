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
#include "engine.h"
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
float supernovae_ii_get_number_per_unit_mass(const struct supernovae_ii *snii,
                                             float m1, float m2) {
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
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 * @param yields The elements ejected (needs to be allocated).
 */
void supernovae_ii_get_yields_from_integral(const struct supernovae_ii *snii,
                                            float log_m1, float log_m2,
                                            float *yields) {

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    float yields_1 = interpolate_1d(&snii->integrated.yields[i], log_m1);
    float yields_2 = interpolate_1d(&snii->integrated.yields[i], log_m2);

    yields[i] = yields_2 - yields_1;
  }
};

/**
 * @brief Get the SNII yields per mass.
 *
 * @param snii The #supernovae_ii model.
 * @param log_m The mass in log.
 * @param yields The elements ejected (needs to be allocated).
 */
void supernovae_ii_get_yields_from_raw(const struct supernovae_ii *snii,
                                       float log_m, float *yields) {

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    yields[i] = interpolate_1d(&snii->raw.yields[i], log_m);
  }
};

/**
 * @brief Get the ejected mass (non processed) per mass unit.
 *
 * @param snii The #supernovae_ii model.
 * @param log_m1 The lower mass in log.
 * @param log_m2 The upper mass in log.
 *
 * @return mass_ejected_processed The mass of non processsed elements.
 */
float supernovae_ii_get_ejected_mass_fraction_non_processed_from_integral(
    const struct supernovae_ii *snii, float log_m1, float log_m2) {

  float mass_ejected_1 =
      interpolate_1d(&snii->integrated.ejected_mass_non_processed, log_m1);
  float mass_ejected_2 =
      interpolate_1d(&snii->integrated.ejected_mass_non_processed, log_m2);

  return mass_ejected_2 - mass_ejected_1;
};

/**
 * @brief Get the ejected mass (non processed) per mass unit.
 *
 * @param snii The #supernovae_ii model.
 * @param log_m The mass in log.
 *
 * @return The mass of non processsed elements.
 */
float supernovae_ii_get_ejected_mass_fraction_non_processed_from_raw(
    const struct supernovae_ii *snii, float log_m) {

  return interpolate_1d(&snii->raw.ejected_mass_non_processed, log_m);
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
float supernovae_ii_get_ejected_mass_fraction_processed_from_integral(
    const struct supernovae_ii *snii, float log_m1, float log_m2) {

  float mass_ejected_1 =
      interpolate_1d(&snii->integrated.ejected_mass_processed, log_m1);
  float mass_ejected_2 =
      interpolate_1d(&snii->integrated.ejected_mass_processed, log_m2);

  return mass_ejected_2 - mass_ejected_1;
};

/**
 * @brief Get the ejected mass (processed) per mass.
 *
 * @param snii The #supernovae_ii model.
 * @param log_m The mass in log.
 *
 * @return mass_ejected The mass of non processsed elements.
 */
float supernovae_ii_get_ejected_mass_fraction_processed_from_raw(
    const struct supernovae_ii *snii, float log_m) {

  return interpolate_1d(&snii->raw.ejected_mass_processed, log_m);
};

/**
 * @brief Get the supernova energy for a given progenitor stellar mass.
 *
 * @param The progenitor mass in units of solar mass.
 *
 * @return the energy released in ergs/1e51.
 */
float supernovae_ii_get_energy_from_progenitor_mass(
    const struct supernovae_ii *snii, float mass) {

  float log_mass = log10(mass);
  return interpolate_1d(&snii->energy_per_progenitor_mass, log_mass);
};

/**
 * @brief Read an array of SNII yields from the table.
 *
 * @param snii The #supernovae_ii model.
 * @param interp_raw Interpolation data to initialize (raw).
 * @param interp_int Interpolation data to initialize (integrated).
 * @param sm * The #stellar_model.
 * @param group_id The HDF5 group id where to read from.
 * @param hdf5_dataset_name The dataset name to read.
 * @param previous_count Number of element in the previous array read.
 * @param interpolation_size Number of element to keep in the interpolation
 * data.
 */
void supernovae_ii_read_yields_array(
    struct supernovae_ii *snii, struct interpolation_1d *interp_raw,
    struct interpolation_1d *interp_int, const struct stellar_model *sm,
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

  /* Initialize the raw interpolation */
  interpolate_1d_init(interp_raw, log10(snii->mass_min), log10(snii->mass_max),
                      interpolation_size, log_mass_min, step_size, count, data,
                      boundary_condition_zero);

  initial_mass_function_integrate(&sm->imf, data, count, log_mass_min,
                                  step_size);
  // TODO: decrease count in order to keep the same distance between points

  /* Initialize the integrated interpolation */
  interpolate_1d_init(interp_int, log10(snii->mass_min), log10(snii->mass_max),
                      interpolation_size, log_mass_min, step_size, count, data,
                      boundary_condition_const);

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
 * @param restart Are we restarting the simulation? (Is params NULL?)
 */
void supernovae_ii_read_yields(struct supernovae_ii *snii,
                               struct swift_params *params,
                               const struct stellar_model *sm,
                               const int restart) {

  hid_t file_id, group_id;

  hsize_t previous_count = 0;

  if (!restart) {
    snii->interpolation_size = parser_get_opt_param_int(
        params, "GEARSupernovaeII:interpolation_size", 200);
  }

  /* Open IMF group */
  h5_open_group(sm->yields_table, "Data/SNII", &file_id, &group_id);

  /* Do all the elements */
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {

    /* Get the element name */
    const char *name = stellar_evolution_get_element_name(sm, i);

    /* Read the array */
    supernovae_ii_read_yields_array(
        snii, &snii->raw.yields[i], &snii->integrated.yields[i], sm, group_id,
        name, &previous_count, snii->interpolation_size);
  }

  /* Read the mass ejected */
  supernovae_ii_read_yields_array(snii, &snii->raw.ejected_mass_processed,
                                  &snii->integrated.ejected_mass_processed, sm,
                                  group_id, "Ej", &previous_count,
                                  snii->interpolation_size);

  /* Read the mass ejected of non processed gas */
  supernovae_ii_read_yields_array(snii, &snii->raw.ejected_mass_non_processed,
                                  &snii->integrated.ejected_mass_non_processed,
                                  sm, group_id, "Ejnp", &previous_count,
                                  snii->interpolation_size);

  /* Cleanup everything */
  h5_close_group(file_id, group_id);
};

/**
 * @brief Reads the supernovae II parameters from tables.
 *
 * @param snii The #supernovae_ii model.
 * @param params The simulation parameters.
 * @param filename The filename of the chemistry table.
 */
void supernovae_ii_read_from_tables(struct supernovae_ii *snii,
                                    struct swift_params *params,
                                    const char *filename) {

  hid_t file_id, group_id;

  /* Open IMF group */
  h5_open_group(filename, "Data/SNII", &file_id, &group_id);

  /* Read the minimal mass of a supernovae */
  io_read_attribute(group_id, "Mmin", FLOAT, &snii->mass_min);

  /* Read the maximal mass of a supernovae */
  io_read_attribute(group_id, "Mmax", FLOAT, &snii->mass_max);

  /* Cleanup everything */
  h5_close_group(file_id, group_id);
}

/**
 * @brief Reads the supernovae II energy parameters from tables.
 *
 * @param snii The #supernovae_ii model.
 * @param params The simulation parameters.
 * @param filename The filename of the chemistry table.
 */
void supernovae_ii_read_energy_from_tables(struct supernovae_ii *snii,
                                           struct swift_params *params,
                                           const char *filename) {

  const char *hdf5_dataset_name = "Energy";
  hid_t file_id, group_id;

  /* Open IMF group */
  h5_open_group(filename, "Data/SNIIEnergy", &file_id, &group_id);

  /* Now let's get the number of elements */
  /* Open attribute */
  const hid_t h_dataset = H5Dopen(group_id, hdf5_dataset_name, H5P_DEFAULT);
  if (h_dataset < 0)
    error("Error while opening attribute '%s'", hdf5_dataset_name);

  /* Get the number of elements */
  hsize_t count = io_get_number_element_in_dataset(h_dataset);

  /* Read the minimal energy */
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
    error("Failed to allocate the SNII energy for %s.", hdf5_dataset_name);

  /* Read the dataset */
  io_read_array_dataset(group_id, hdf5_dataset_name, FLOAT, data, count);

  /* Initialize the raw interpolation */
  interpolate_1d_init(&snii->energy_per_progenitor_mass, log10(snii->mass_min),
                      log10(snii->mass_max), snii->interpolation_size,
                      log_mass_min, step_size, count, data,
                      boundary_condition_zero);

  /* Cleanup the memory */
  free(data);

  /* Cleanup everything */
  h5_close_group(file_id, group_id);
}

/**
 * @brief Initialize the #supernovae_ii structure.
 *
 * @param snii The #supernovae_ii model.
 * @param params The simulation parameters.
 * @param sm The #stellar_model.
 * @param us The unit system.
 */
void supernovae_ii_init(struct supernovae_ii *snii, struct swift_params *params,
                        const struct stellar_model *sm,
                        const struct unit_system *us) {

  /* Read the parameters from the tables */
  supernovae_ii_read_from_tables(snii, params, sm->yields_table);

  /* Read the supernovae yields */
  supernovae_ii_read_yields(snii, params, sm, /* restart */ 0);

  /* Get the IMF parameters */
  snii->exponent = initial_mass_function_get_exponent(&sm->imf, snii->mass_min,
                                                      snii->mass_max);
  snii->coef_exp = initial_mass_function_get_coefficient(
      &sm->imf, snii->mass_min, snii->mass_max);
  snii->coef_exp /= snii->exponent;

  /* Read the energy parameters from the tables */
  supernovae_ii_read_energy_from_tables(snii, params, sm->yields_table);

  /* Supernovae energy */
  double e_feedback =
      parser_get_param_double(params, "GEARFeedback:supernovae_energy_erg");
  e_feedback /= units_cgs_conversion_factor(us, UNIT_CONV_ENERGY);
  snii->energy_per_supernovae = e_feedback;
}

/**
 * @brief Write a supernovae_ii struct to the given FILE as a stream of bytes.
 *
 * Here we are only writing the arrays, everything else has been copied in the
 * feedback.
 *
 * @param snii the struct
 * @param stream the file stream
 * @param sm The #stellar_model.
 */
void supernovae_ii_dump(const struct supernovae_ii *snii, FILE *stream,
                        const struct stellar_model *sm) {}

/**
 * @brief Restore a supernovae_ii struct from the given FILE as a stream of
 * bytes.
 *
 * Here we are only writing the arrays, everything else has been copied in the
 * feedback.
 *
 * @param snii the struct
 * @param stream the file stream
 * @param sm The #stellar_model.
 */
void supernovae_ii_restore(struct supernovae_ii *snii, FILE *stream,
                           const struct stellar_model *sm) {

  /* Read the supernovae yields (and apply the units) */
  supernovae_ii_read_yields(snii, NULL, sm, /* restart */ 1);

  /* Get the IMF parameters */
  snii->exponent = initial_mass_function_get_exponent(&sm->imf, snii->mass_min,
                                                      snii->mass_max);
  snii->coef_exp = initial_mass_function_get_coefficient(
      &sm->imf, snii->mass_min, snii->mass_max);
  snii->coef_exp /= snii->exponent;

  /* Read the energy parameters from the tables */
  supernovae_ii_read_energy_from_tables(snii, NULL, sm->yields_table);
}

/**
 * @brief Clean the allocated memory.
 *
 * @param snii the #supernovae_ii.
 */
void supernovae_ii_clean(struct supernovae_ii *snii) {

  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    interpolate_1d_free(&snii->integrated.yields[i]);
    interpolate_1d_free(&snii->raw.yields[i]);
  }

  interpolate_1d_free(&snii->integrated.ejected_mass_processed);
  interpolate_1d_free(&snii->raw.ejected_mass_processed);
  interpolate_1d_free(&snii->integrated.ejected_mass_non_processed);
  interpolate_1d_free(&snii->raw.ejected_mass_non_processed);
}
