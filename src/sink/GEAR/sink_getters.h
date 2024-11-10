/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2024 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
 *******************************************************************************/
#ifndef SWIFT_GEAR_SINK_GETTERS_H
#define SWIFT_GEAR_SINK_GETTERS_H

#include "cosmology.h"
#include "sink_part.h"

/**
 * @file src/sink/GEAR/sink_getter.h
 * @brief Getter functions for GEAR sink scheme to avoid exposing
 * implementation details to the outer world. Keep the code clean and lean.
 */

/**
 * @brief Get the sink age in interal units.
 *
 * @param sink The #sink.
 * @param with_cosmology If we run with cosmology.
 * @param cosmo The co Birth scale-factor of the star.
 * @param with_cosmology If we run with cosmology.
 */

__attribute__((always_inline)) INLINE double
sink_get_sink_age(struct sink* restrict sink, const int with_cosmology,
		  const struct cosmology* cosmo, const int time) {
  double sink_age;
  if (with_cosmology) {

    /* Deal with rounding issues */
    if (sink->birth_scale_factor >= cosmo->a) {
      sink_age = 0.;
    } else {
      sink_age = cosmology_get_delta_time_from_scale_factors(cosmo, sink->birth_scale_factor, cosmo->a);
    }
  } else {
    sink_age = time - sink->birth_time;
  }
  return sink_age;
}

#endif /* SWIFT_GEAR_SINK_GETTERS_H */
