/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2025 Zachary Rey (zachary.rey@epfl.ch)
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

/* Include Header*/
#include "stellar_wind.h"

#include "interpolation.h"

#include <math.h>
#include <stdio.h>
// #include "unit.h"

// TODO: Do we want to print properties of the stellar wind?

/**
 * @brief Read an array of stellar wind from the table.
 *
 * @param sw The #stellar_wind model.
 * @param interp Interpolation data to initialize.
 * @param sm * The #stellar_model.
 * @param group_id The HDF5 group id where to read from.
 * @param hdf5_dataset_name The dataset name to read.
 * @param previous_count Number of element in the previous array read.
 * @param interpolation_size_m Number of element to keep in the mass
 * interpolation data.
 * @param interpolation_size_z Number of element to keep in the metallicity
 * interpolation data.
 */
void stellar_wind_read_yields_array(
    struct stellar_wind *sw, struct interpolation_2d *interp,
    const struct stellar_model *sm, hid_t group_id,
    const char *hdf5_dataset_name, hsize_t *previous_count,
    int interpolation_size_m, int interpolation_size_z) {

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
  io_read_attribute(h_dataset, "m0", FLOAT, &log_mass_min);

  /* Read the step size in mass array (log step) */
  float step_size_m = 0;
  io_read_attribute(h_dataset, "dm", FLOAT, &step_size_m);

  /* Read the number of element in mass array */
  int Nm = 0;
  io_read_attribute(h_dataset, "nm", INT, &Nm);

  /* Read the minimal metallicity (in log) */
  float log_metallicity_min = 0;
  io_read_attribute(h_dataset, "z0", FLOAT, &log_metallicity_min);

  /* Read the step size in metallicity array (log step) */
  float step_size_z = 0;
  io_read_attribute(h_dataset, "dz", FLOAT, &step_size_z);

  /* Read the number of element in metallicity array */
  int Nz = 0;
  io_read_attribute(h_dataset, "nz", INT, &Nz);

  /* Close the attribute */
  H5Dclose(h_dataset);

  /* Allocate the memory */
  double *data = (double *)malloc(sizeof(double) * Nz * Nm);
  if (data == NULL)
    error("Failed to allocate the SW yields for %s.", hdf5_dataset_name);

  /* Read the dataset */
  io_read_array_dataset(group_id, hdf5_dataset_name, DOUBLE, data, count);

  /* Initialize the raw interpolation */
  float log_mass_max = log_mass_min + step_size_m * Nm;
  float log_metallicity_max = log_metallicity_min + step_size_z * Nz;

  interpolate_2d_init(interp, log_metallicity_min, log_metallicity_max,
                      interpolation_size_z, log_mass_min, log_mass_max,
                      interpolation_size_m, log_metallicity_min, log_mass_min,
                      step_size_z, step_size_m, Nz, Nm, data,
                      boundary_condition_const);

  /* Cleanup the memory */
  free(data);
}

/**
 * @brief Read the SW yields from the table.
 *
 * The tables are in [erg/yr] units at the end of this function. TODO: convert
 * into internal units
 *
 * @param sw The #stellar_wind model.
 * @param params The simulation parameters.
 * @param sm The #stellar_model.
 * @param restart Are we restarting the simulation? (Is params NULL?)
 */
void stellar_wind_read_yields(struct stellar_wind *sw,
                              struct swift_params *params,
                              const struct stellar_model *sm,
                              const int restart) {

  hid_t file_id, group_id;

  hsize_t previous_count = 0;

  if (!restart) {
    sw->interpolation_size_m = parser_get_opt_param_int(
        params, "GEARStellar_wind:interpolation_size_mass", 200);
    sw->interpolation_size_z = parser_get_opt_param_int(
        params, "GEARStellar_wind:interpolation_size_metallicity", 110);
  }

  /* Open IMF group */
  h5_open_group(sm->yields_table, "Data/SW/MetallicityDependent", &file_id,
                &group_id);

  /* Read the energy*/
  stellar_wind_read_yields_array(
      sw, &sw->raw.ejected_energy, sm, group_id, "Energy", &previous_count,
      sw->interpolation_size_m, sw->interpolation_size_z);

  /* Read the integrated energy*/
  stellar_wind_read_yields_array(
      sw, &sw->integrated.ejected_energy_per_progenitor_mass, sm, group_id,
      "Integrated_Energy", &previous_count, sw->interpolation_size_m,
      sw->interpolation_size_z);

  /* Read the mass-loss*/
  stellar_wind_read_yields_array(
      sw, &sw->raw.mass_loss, sm, group_id, "Mass_Loss", &previous_count,
      sw->interpolation_size_m, sw->interpolation_size_z);

  /* Read the integrated mass-loss*/
  stellar_wind_read_yields_array(
      sw, &sw->integrated.mass_loss_per_progenitor_mass, sm, group_id,
      "Integrated_Mass_Loss", &previous_count, sw->interpolation_size_m,
      sw->interpolation_size_z);

  /* Cleanup everything */
  h5_close_group(file_id, group_id);
};

/**
 * @brief Initialize the #stellar_wind structure.
 *
 * @param sw The #stellar_wind model.
 * @param params The simulation parameters.
 * @param sm The #stellar_model.
 * @param us The unit system.
 */
void stellar_wind_init(struct stellar_wind *sw, struct swift_params *params,
                       const struct stellar_model *sm,
                       const struct unit_system *us) {

  /* Read the stellar wind yields */
  stellar_wind_read_yields(sw, params, sm, 0);
}

/**
 * @brief Restore a stellar_wind struct from the given FILE as a stream of
 * bytes.
 *
 * @param sw The #stellar_wind model.
 * @param params The simulation parameters.
 * @param sm The #stellar_model.
 */
void stellar_wind_restore(struct stellar_wind *sw, FILE *stream,
                          const struct stellar_model *sm) {

  stellar_wind_read_yields(sw, NULL, sm, 1);
}

/**
 * @brief Get the ejected energy given a discrete mass.
 *
 * @param sw The #stellar_wind model.
 * @param log_m The upper mass in log.
 * @param log_z The metallicity in log.
 *
 * @return energy per unit time in [erg/yr].
 */
double stellar_wind_get_ejected_energy(const struct stellar_wind *sw,
                                       float log_m, float log_z) {
  return pow(10, interpolate_2d(&sw->raw.ejected_energy, log_z, log_m));
};

/**
 * @brief Get the ejected energy per progenitor mass.
 *
 * @param sw The #stellar_wind model.
 * @param log_m The upper mass in log.
 * @param log_z The metallicity in log.
 *
 * @return energy per progenitor mass per unit time in [erg/yr].
 */
double stellar_wind_get_ejected_energy_IMF(const struct stellar_wind *sw,
                                           float log_m, float log_z) {
  return pow(10,
             interpolate_2d(&sw->integrated.ejected_energy_per_progenitor_mass,
                            log_z, log_m));
};

/**
 * @brief Get the ejected mass given a discrete mass.
 *
 * @param sw The #stellar_wind model.
 * @param log_m The upper mass in log.
 * @param log_z The metallicity in log.
 *
 * @return mass per unit time in [Msol/yr].
 */
double stellar_wind_get_ejected_mass(const struct stellar_wind *sw, float log_m,
                                     float log_z) {
  return pow(10, interpolate_2d(&sw->raw.mass_loss, log_z, log_m));
};

/**
 * @brief Get the ejected mass per progenitor mass.
 *
 * @param sw The #stellar_wind model.
 * @param log_m The upper mass in log.
 * @param log_z The metallicity in log.
 *
 * @return mass per progenitor mass per unit time in [Msol/yr].
 */
double stellar_wind_get_ejected_mass_IMF(const struct stellar_wind *sw,
                                         float log_m, float log_z) {
  return pow(10, interpolate_2d(&sw->integrated.mass_loss_per_progenitor_mass,
                                log_z, log_m));
};

/**
 * @brief Clean the allocated memory.
 *
 * @param sw the #stellar_wind.
 */
void stellar_wind_clean(struct stellar_wind *sw) {

  interpolate_2d_free(&sw->integrated.ejected_energy_per_progenitor_mass);
  interpolate_2d_free(&sw->raw.ejected_energy);
  interpolate_2d_free(&sw->raw.mass_loss);
  interpolate_2d_free(&sw->integrated.mass_loss_per_progenitor_mass);
}
