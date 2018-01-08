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
#ifndef SWIFT_COOLING_STRUCT_NONE_H
#define SWIFT_COOLING_STRUCT_NONE_H

/* include grackle */
#include <grackle.h>

/**
 * @file src/cooling/none/cooling_struct.h
 * @brief Empty infrastructure for the cases without cooling function
 */

/**
 * @brief Properties of the cooling function.
 */
struct cooling_function_data {

  /* Filename of the Cloudy Table */
  char cloudy_table[200];

  /* Enable/Disable UV backgroud */
  int uv_background;

  /* Redshift to use for the UV backgroud (-1 to use cosmological one) */
  double redshift;

  /* Density Threshold for the shielding */
  double density_self_shielding;

  /* unit system */
  code_units units;

  /* grackle chemistry data */
  chemistry_data chemistry;
};

/**
 * @brief Properties of the cooling stored in the particle data
 */
struct cooling_xpart_data {

  /*! Energy radiated away by this particle since the start of the run */
  float radiated_energy;
};

#endif /* SWIFT_COOLING_STRUCT_NONE_H */
