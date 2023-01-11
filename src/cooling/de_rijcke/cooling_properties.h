//
// Created by yuyttenh on 10/01/23.
//

#ifndef SWIFTSIM_COOLING_PROPERTIES_DE_RIJCKE_H
#define SWIFTSIM_COOLING_PROPERTIES_DE_RIJCKE_H

#define de_rijcke_table_path_name_length 500

/**
 * @brief struct containing cooling tables
 */
struct cooling_tables {
  /* array of temperatures */
  float *temperature;

  /* array of cooling rates */
  float *cooling_rate;
};

/**
 * @brief Properties of the cooling function.
 */
struct cooling_function_data {

  /*! Cooling tables */
  struct cooling_tables table;

  /*! Filepath to the directory containing the HDF5 cooling tables */
  char cooling_table_path[de_rijcke_table_path_name_length];

  /*! Conversion factor from internal units to cgs for density */
  double conv_factor_density_to_cgs;

  /*! Conversion factor from internal units from cgs for internal energy
   * derivative */
  double conv_factor_energy_rate_from_cgs;

  /*! Inverse of the proton mass in cgs units [g^-1] */
  double proton_mass_cgs_inv;

  /*! Minimally allowed internal energy of the particles */
  float min_energy;

  /*! Use rapid cooling? */
  int rapid_cooling;
};

#endif  // SWIFTSIM_COOLING_PROPERTIES_DE_RIJCKE_H
