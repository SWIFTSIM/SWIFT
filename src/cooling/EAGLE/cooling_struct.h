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
#include <stdbool.h>
#include <time.h>

// EAGLE defined constants
#define eagle_element_name_length 20
#define eagle_cmb_temperature 2.728
#define eagle_compton_rate 1.0178e-37 * 2.728 * 2.728 * 2.728 * 2.728
#define eagle_metal_cooling_on 1
#define eagle_max_iterations 15
#define eagle_proton_mass_cgs 1.6726e-24

/*
 * @brief struct containing radiative and heating rates
 */
// struct radiative_rates{
//  float Cooling[NUMBER_OF_CTYPES];
//  float Heating[NUMBER_OF_HTYPES];
//  float TotalCoolingRate;
//  float TotalHeatingRate;
//  float CoolingTime; /*Cooling time in Myr*/
//  float HeatingTime; /*Heating time in Myr*/
//};

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
  struct cooling_tables_redshift_invariant photoionisation_cooling;
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

  int N_Redshifts;
  int N_nH;
  int N_Temp;
  int N_He;
  int N_Elements;
  int N_SolarAbundances;

  float *Redshifts;
  float *nH;
  float *Temp;
  float *HeFrac;
  float *Therm;
  float *SolarAbundances;
  float *SolarElectronAbundance;
  float *ElementAbundance_SOLARM1;
  double *solar_abundances;
  char **ElementNames;
  char **SolarAbundanceNames;
  int *ElementNamePointers;
  int *SolarAbundanceNamePointers;

  /*! Constant multiplication factor for time-step criterion */
  float cooling_tstep_mult;
  float Redshift;
  float min_energy;
  float cooling_rate;
  double internal_energy_scale;
  double number_density_scale;
  double temperature_scale;
  double power_scale;
  char cooling_table_path[500];
  float reionisation_redshift;
  float calcium_over_silicon_ratio;
  float sulphur_over_silicon_ratio;

  double delta_u;
};

/**
 * @brief Properties of the cooling stored in the extended particle data.
 */
struct cooling_xpart_data {
  float radiated_energy;
};

#endif /* SWIFT_COOLING_STRUCT_EAGLE_H */
