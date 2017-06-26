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
#ifndef SWIFT_HYDRO_PROPERTIES
#define SWIFT_HYDRO_PROPERTIES

/**
 * @file hydro_properties.h
 * @brief Contains all the constants and parameters of the hydro scheme
 */

/* Config parameters. */
#include "../config.h"

#if defined(HAVE_HDF5)
#include <hdf5.h>
#endif

/* Local includes. */
#include "parser.h"

/**
 * @brief Contains all the constants and parameters of the hydro scheme
 */
struct hydro_props {

  /*! Resolution parameter */
  float eta_neighbours;

  /*! Target weightd number of neighbours (for info only)*/
  float target_neighbours;

  /*! Smoothing length tolerance */
  float h_tolerance;

  /*! Tolerance on neighbour number  (for info only)*/
  float delta_neighbours;

  /*! Maximal smoothing length */
  float h_max;

  /*! Maximal number of iterations to converge h */
  int max_smoothing_iterations;

  /*! Time integration properties */
  float CFL_condition;

  /*! Maximal change of h over one time-step */
  float log_max_h_change;
};

void hydro_props_print(const struct hydro_props *p);
void hydro_props_init(struct hydro_props *p, const struct swift_params *params);

#if defined(HAVE_HDF5)
void hydro_props_print_snapshot(hid_t h_grpsph, const struct hydro_props *p);
#endif

#endif /* SWIFT_HYDRO_PROPERTIES */
