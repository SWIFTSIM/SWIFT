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
#ifndef SWIFT_COOLING_PROPERTIES_GRACKLE_H
#define SWIFT_COOLING_PROPERTIES_GRACKLE_H

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

  /*! Redshift to use for the UV backgroud (-1 to use cosmological one) */
  double redshift;

  /*! unit system */
  code_units units;

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
