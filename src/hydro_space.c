/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2017 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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

#include "hydro_space.h"

#include "space.h"

/**
 * @brief Initialize the extra space information needed for some hydro schemes.
 *
 * @param hs #hydro_space to initialize.
 * @param s #space containing the hydro space.
 */
#ifdef SHADOWSWIFT
void hydro_space_init(struct hydro_space *hs, const struct space *s,
                      struct swift_params *params) {
#if (SHADOWSWIFT_BC == INFLOW_BC || SHADOWSWIFT_BC == RADIAL_INFLOW_BC)
  if (!s->periodic) {
    hs->density =
        parser_get_param_float(params, "InitialConditions:inflow_density");
    hs->velocity =
        parser_get_param_float(params, "InitialConditions:inflow_velocity");
    hs->pressure =
        parser_get_param_float(params, "InitialConditions:inflow_pressure");
  }
#else
  hs->density = 0.f;
  hs->velocity = 0.f;
  hs->pressure = 0.f;
#endif
  hs->center[0] = 0.5 * s->dim[0];
  hs->center[1] = 0.5 * s->dim[1];
  hs->center[2] = 0.5 * s->dim[2];
}
#else
void hydro_space_init(struct hydro_space* hs, const struct space* s,
                      struct swift_params* params) {}
#endif
