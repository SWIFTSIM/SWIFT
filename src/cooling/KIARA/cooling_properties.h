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
#ifndef SWIFT_COOLING_PROPERTIES_KIARA_H
#define SWIFT_COOLING_PROPERTIES_KIARA_H

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

  /*! Redshift to use for the UV backgroud (-1 to use cosmological one) */
  double redshift;

  /*! unit system */
  code_units units;

  /*! k_Boltz/m_p plus conversion factor for converting u<->T */
  double temp_to_u_factor;

  /*! conversion unit factor for rate of change of thermal energy */
  double dudt_units;

  /*! grackle chemistry data */
  chemistry_data chemistry;

  /*! Enable/Disable metal cooling */
  int with_metal_cooling;

  /*! User provide volumetric heating rates */
  int provide_volumetric_heating_rates;

  /*! User provide specific heating rates */
  int provide_specific_heating_rates;

  /*! Self shielding method (1 -> 3 for grackle's ones, 0 for none and -1 for
   * GEAR) */
  int self_shielding_method;

  /*! convergence limit for first init */
  float convergence_limit;

  /*! number of step max for first init */
  int max_step;

  /*! Duration for switching off cooling after an event (e.g. supernovae) */
  double thermal_time;

  /*! track dust growth and destruction (only available in KIARA) */
  int use_grackle_dust_evol;

  /*! track H2 formation; this is set within the code based on selection options */
  int use_grackle_h2_form;

  /*! G0 conversion factor, scales to MW value based on galaxy-wide sSFR */
  double G0_factor;

  /*! Dust parameters; see sample yml file */
  double dust_destruction_eff;
  double dust_sne_coeff;
  double dust_sne_shockspeed;
  double dust_grainsize;
  double dust_growth_densref;
  double dust_growth_tauref;

  /*! For subgrid model (eg KIARA) need a subgrid ISM fraction */
  double cold_ISM_frac;
};

#endif /* SWIFT_COOLING_PROPERTIES_KIARA_H */
