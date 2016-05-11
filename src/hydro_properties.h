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

/* Config parameters. */
#include "../config.h"

/* Local includes. */
#include "const.h"
#include "parser.h"

/**
 * @brief Contains all the constants and parameters of the hydro scheme
 */
struct hydro_props {

  /* Kernel properties */
  float eta_neighbours;
  float target_neighbours;
  float delta_neighbours;

  /* Kernel properties */
  int max_smoothing_iterations;

  /* Time integration properties */
  float CFL_condition;
  float log_max_h_change;

/* Viscosity parameters */
#ifdef GADGET_SPH
  float const_viscosity_alpha;
#endif
};

void hydro_props_print(const struct hydro_props *p);
void hydro_props_init(struct hydro_props *p, const struct swift_params *params);

#endif /* SWIFT_HYDRO_PROPERTIES */
