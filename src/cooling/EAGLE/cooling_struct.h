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

// EAGLE defined constants
#define eagle_element_name_length 20
#define eagle_metal_cooling_on 1
#define eagle_max_iterations 15

/*
 * @brief struct containing cooling tables independent of redshift
 *
 */
struct cooling_tables_redshift_invariant {
  float *metal_heating;      // size: [1][cooling_N_Temp][cooling_N_nH];
  float *H_plus_He_heating;  // size:
                             // [1][cooling_N_He][cooling_N_Temp][cooling_N_nH];
  float *
      H_plus_He_electron_abundance;  // size:
                                     // [1][cooling_N_He][cooling_N_Temp][cooling_N_nH];
  float *temperature;  // size: [1][cooling_N_He][cooling_N_Temp][cooling_N_nH];
  float *electron_abundance;  // size: [1][cooling_N_Temp][cooling_N_nH];
};

/*
 * @brief struct containing cooling tables depending on redshift
 *
 */
struct cooling_tables {
  float *metal_heating;  // size:
                         // [cooling_N_Redshifts][cooling_N_Temp][cooling_N_nH];
  float *
      H_plus_He_heating;  // size:
                          // [cooling_N_Redshifts][cooling_N_He][cooling_N_Temp][cooling_N_nH];
  float *
      H_plus_He_electron_abundance;  // size:
                                     // [cooling_N_Redshifts][cooling_N_He][cooling_N_Temp][cooling_N_nH];
  float *
      temperature;  // size:
                    // [cooling_N_Redshifts][cooling_N_He][cooling_N_Temp][cooling_N_nH];
  float *
      electron_abundance;  // size:
                           // [cooling_N_Redshifts][cooling_N_Temp][cooling_N_nH];
};

/*
 * @brief struct containing structs of various eagle cooling tables
 *
 */
struct eagle_cooling_table {
  struct cooling_tables_redshift_invariant no_compton_cooling;
  struct cooling_tables_redshift_invariant photodissociation_cooling;
  struct cooling_tables_redshift_invariant collisional_cooling;
  struct cooling_tables element_cooling;
};

/**
 * @brief Properties of the cooling function.
 */
struct cooling_function_data {

  /* Cooling table variables */
  struct eagle_cooling_table table;

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


  /* Arrays containing names of tracked elements. Used for 
   * reading in the contributions to the cooling rate from
   * each element */
  char **ElementNames;
  char **SolarAbundanceNames;

  /* Normalisation constants that are frequently used */
  double internal_energy_scale;
  double number_density_scale;
  double temperature_scale;
  double power_scale;

  /* filepath to EAGLE cooling tables */
  char cooling_table_path[500];

  /* Some constants read in from yml file relevant to EAGLE cooling */
  float reionisation_redshift;
  float calcium_over_silicon_ratio;
  float sulphur_over_silicon_ratio;

  /* Helium reionisation parameters */
  int he_reion_flag;
  float he_reion_ev_pH;
  float he_reion_z_center;
  float he_reion_z_sigma;

  /* Proton mass in cgs */
  double proton_mass_cgs;
  double T_CMB_0;
  double compton_rate_cgs;

};

/**
 * @brief Properties of the cooling stored in the extended particle data.
 */
struct cooling_xpart_data {
  float radiated_energy;
};

#endif /* SWIFT_COOLING_STRUCT_EAGLE_H */
