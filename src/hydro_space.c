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

#include "error.h"
#include "space.h"

/**
 * @brief Initialize the extra space information needed for some hydro schemes.
 *
 * @param hs #hydro_space to initialize.
 * @param s #space containing the hydro space.
 */
#ifdef SHADOWFAX_SPH
__attribute__((always_inline)) INLINE void hydro_space_init(
    struct hydro_space *hs, const struct space *s) {

  if (s->periodic) {
    hs->anchor[0] = -0.5f * s->dim[0];
    hs->anchor[1] = -0.5f * s->dim[1];
    hs->anchor[2] = -0.5f * s->dim[2];
    hs->side[0] = 2.0f * s->dim[0];
    hs->side[1] = 2.0f * s->dim[1];
    hs->side[2] = 2.0f * s->dim[2];
  } else {
    hs->anchor[0] = 0.0f;
    hs->anchor[1] = 0.0f;
    hs->anchor[2] = 0.0f;
    hs->side[0] = s->dim[0];
    hs->side[1] = s->dim[1];
    hs->side[2] = s->dim[2];
  }

  message(
      "Initialised hydro space with anchor [%g, %g, %g] and side [%g, %g, %g].",
      hs->anchor[0], hs->anchor[1], hs->anchor[2], hs->side[0], hs->side[1],
      hs->side[2]);
}
#else
void hydro_space_init(struct hydro_space *hs, const struct space *s) {}
#endif
