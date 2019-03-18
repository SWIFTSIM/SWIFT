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

/**
 * @brief struct containing cooling tables
 */
struct cooling_tables {

  /* array of heating rates due to metals */
  float *metal_heating;

  /* array of heating rates due to hydrogen and helium */
  float *H_plus_He_heating;

  /* array of electron abundances due to hydrogen and helium */
  float *H_plus_He_electron_abundance;

  /* array of temperatures */
  float *temperature;

  /* array of electron abundances due to metals */
  float *electron_abundance;
};

/**
 * @brief Properties of the cooling function.
 */
struct cooling_function_data {

  /*! Cooling tables */
  struct cooling_tables table;

  /*! Redshift bins */
  float *Redshifts;

  /*! Hydrogen number density bins */
  float *nH;

  /*! Temperature bins */
  float *Temp;

  /*! Helium fraction bins */
  float *HeFrac;

  /*! Internal energy bins */
  float *Therm;

  /*! Mass fractions of elements for solar abundances (from the tables) */
  float *SolarAbundances;

  /*! Inverse of the solar mass fractions */
  float *SolarAbundances_inv;

  /*! Filepath to the directory containing the HDF5 cooling tables */
  char cooling_table_path[eagle_table_path_name_length];

  /*! Redshift of H reionization */
  float H_reion_z;

  /*! H reionization energy in CGS units */
  float H_reion_heat_cgs;

  /*! Have we already done H reioisation? */
  int H_reion_done;

  /*! Ca over Si abundance divided by the solar ratio for these elements */
  float Ca_over_Si_ratio_in_solar;

  /*! S over Si abundance divided by the solar ratio for these elements */
  float S_over_Si_ratio_in_solar;

  /*! Redshift of He reionization */
  float He_reion_z_centre;

  /*! Spread of the He reionization */
  float He_reion_z_sigma;

  /*! He reionization energy in CGS units */
  float He_reion_heat_cgs;

  /*! Internal energy conversion from internal units to CGS (for quick access)
   */
  double internal_energy_to_cgs;

  /*! Internal energy conversion from CGS to internal units (for quick access)
   */
  double internal_energy_from_cgs;

  /*! Number density conversion from internal units to CGS (for quick access) */
  double number_density_to_cgs;

  /*! Inverse of proton mass in cgs (for quick access) */
  double inv_proton_mass_cgs;

  /*! Temperatur of the CMB at present day (for quick access) */
  double T_CMB_0;

  /*! Compton rate in cgs units */
  double compton_rate_cgs;

  /*! Index of the current redshift along the redshift index of the tables */
  int z_index;

  /*! Distance between the current redshift and table[z_index] */
  float dz;

  /*! Index of the previous tables along the redshift index of the tables */
  int previous_z_index;

  /*! Are we doing Newton-Raphson iterations? */
  int newton_flag;
};

/**
 * @brief Properties of the cooling stored in the extended particle data.
 */
struct cooling_xpart_data {

  /*! Cumulative energy radiated by the particle */
  float radiated_energy;
};

#endif /* SWIFT_COOLING_STRUCT_EAGLE_H */
