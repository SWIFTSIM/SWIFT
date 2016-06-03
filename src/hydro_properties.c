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

/* This object's header. */
#include "hydro_properties.h"

/* Standard headers */
#include <float.h>
#include <math.h>

/* Local headers. */
#include "error.h"
#include "kernel_hydro.h"
#include "hydro.h"

void hydro_props_init(struct hydro_props *p,
                      const struct swift_params *params) {

  /* Kernel properties */
  p->eta_neighbours = parser_get_param_float(params, "SPH:resolution_eta");
  const float eta3 = p->eta_neighbours * p->eta_neighbours * p->eta_neighbours;
  p->target_neighbours = 4.0 * M_PI * kernel_gamma3 * eta3 / 3.0;
  p->delta_neighbours = parser_get_param_float(params, "SPH:delta_neighbours");

  /* Ghost stuff */
  p->max_smoothing_iterations =
      parser_get_param_int(params, "SPH:max_ghost_iterations");

  /* Time integration properties */
  p->CFL_condition = parser_get_param_float(params, "SPH:CFL_condition");
  const float max_volume_change =
      parser_get_param_float(params, "SPH:max_volume_change");
  p->log_max_h_change = logf(powf(max_volume_change, 0.33333333333f));
}

void hydro_props_print(const struct hydro_props *p) {

  message("Hydrodynamic scheme: %s.", SPH_IMPLEMENTATION);
  message("Hydrodynamic kernel: %s with %.2f +/- %.2f neighbours (eta=%f).",
          kernel_name, p->target_neighbours, p->delta_neighbours,
          p->eta_neighbours);
  message("Hydrodynamic integration: CFL parameter: %.4f.", p->CFL_condition);
  message(
      "Hydrodynamic integration: Max change of volume: %.2f "
      "(max|dlog(h)/dt|=%f).",
      powf(expf(p->log_max_h_change), 3.f), p->log_max_h_change);
}
