/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_COOLING_STRUCT_EAGLE_H
#define SWIFT_COOLING_STRUCT_EAGLE_H

#define eagle_table_path_name_length 500

/*
 * @brief struct containing cooling tables
 *
 */
struct cooling_tables {
  // array of heating rates due to metals
  float *metal_heating;

  // array of heating rates due to hydrogen and helium
  float *H_plus_He_heating;

  // array of electron abundances due to hydrogen and helium
  float *H_plus_He_electron_abundance;

  // array of temperatures
  float *temperature;

  // array of electron abundances due to metals
  float *electron_abundance;
};

/**
 * @brief Properties of the cooling function.
 */
struct cooling_function_data {

  /* Cooling table */
  struct cooling_tables table;

  /* Size of table dimensions */
  int N_Redshifts;
  int N_nH;
  int N_Temp;
  int N_He;

  /* Number of metals and solar abundances tracked in EAGLE */
  int N_Elements;
  int N_SolarAbundances;

  /* Arrays of grid values in tables */
  float *Redshifts;
  float *nH;
  float *Temp;
  float *HeFrac;
  float *Therm;

  /* Array of values of solar metal and electron abundance */
  float *SolarAbundances;
  float *SolarElectronAbundance;

  /* Normalisation constants that are frequently used 
   * Multiply by these values to go from internal to cgs
   * units for relevant quantity */
  double internal_energy_scale;
  double number_density_scale;

  /* filepath to EAGLE cooling tables */
  char cooling_table_path[eagle_table_path_name_length];

  /* Some constants read in from yml file relevant to EAGLE cooling */
  float reionisation_redshift;
  float calcium_over_silicon_ratio;
  float sulphur_over_silicon_ratio;

  /* Hydrogen reionisation parameters */
  float h_reion_z_center;

  /* Helium reionisation parameters */
  float he_reion_ev_pH;
  float he_reion_z_center;
  float he_reion_z_sigma;

  /* Proton mass in cgs */
  double proton_mass_cgs;
  double T_CMB_0;
  double compton_rate_cgs;

  /* redshift table indices and offset */
  int z_index, previous_z_index;
  float dz;
  int low_z_index, high_z_index;
};

/**
 * @brief Properties of the cooling stored in the extended particle data.
 */
struct cooling_xpart_data {
  float radiated_energy;
};

#endif /* SWIFT_COOLING_STRUCT_EAGLE_H */
