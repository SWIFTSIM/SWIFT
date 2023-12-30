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
#ifndef SWIFT_COOLING_PROPERTIES_GRACKLE_H
#define SWIFT_COOLING_PROPERTIES_GRACKLE_H

/* skip deprecation warnings. I cleaned old API calls. */
#define OMIT_LEGACY_INTERNAL_GRACKLE_FUNC

/* include grackle */
#include <grackle.h>

/**
 * @file src/cooling/grackle/cooling_properties.h
 * @brief Empty infrastructure for the cases without cooling function
 */

/**
 * @brief Properties of the cooling function.
 */
struct cooling_function_data {

  /*! Filename of the Cloudy Table */
  char cloudy_table[200];

  /*! Enable/Disable UV backgroud */
  int with_uv_background;

  /*! Chemistry network */
  int primordial_chemistry;

  /*! Set the three-body reaction rate (see grackle documentation) */
  int H2_three_body_rate;

  /*! Enable/disable H2 collision-induced emission cooling from Ripamonti & Abel
   * (2004) */
  int H2_cie_cooling;

  /*! Flag to enable H2 formation on dust grains */
  int H2_on_dust;

  /*! The ratio of total dust mass to gas mass in the local Universe. */
  double local_dust_to_gas_ratio;

  /*! Enable/disable CMB temperature floor */
  int cmb_temperature_floor;

  /*! Redshift to use for the UV backgroud (-1 to use cosmological one) */
  double redshift;

  /*! unit system */
  code_units units;

  /*! grackle chemistry data */
  chemistry_data chemistry_data;

  /*! grackle chemistry data storage
   * (needed for local function calls) */
  chemistry_data_storage chemistry_rates;

  /*! Enable/Disable metal cooling */
  int with_metal_cooling;

  /*! Arrays of ionization and heating rates are provided */
  int use_radiative_transfer;

  /*! Grackle RT_heating_rate (in IU) */
  float RT_heating_rate;

  /*! Grackle RT_HI_ionization_rate (in IU) */
  float RT_HI_ionization_rate;

  /*! Grackle RT_HeI_ionization_rate (in IU) */
  float RT_HeI_ionization_rate;

  /*! Grackle RT_HeII_ionization_rate (in IU) */
  float RT_HeII_ionization_rate;

  /*! Grackle RT_H2_dissociation_rate (in IU) */
  float RT_H2_dissociation_rate;

  /*! Volumetric heating rates */
  float volumetric_heating_rates;

  /*! Specific heating rates */
  float specific_heating_rates;

  /*! Hydrogen fraction by mass */
  float HydrogenFractionByMass;

  /*! Self shielding method (1 -> 3 for grackle's ones, 0 for none and -1 for
   * GEAR) */
  int self_shielding_method;

  /*! Self shielding threshold */
  float self_shielding_threshold;

  /*! convergence limit for first init */
  float convergence_limit;

  /*! number of step max for first init */
  int max_step;

  /*! over relaxation parameter */
  float omega;

  /*! Duration for switching off cooling after an event (e.g. supernovae) */
  double thermal_time;
};

#endif /* SWIFT_COOLING_PROPERTIES_GRACKLE_H */
