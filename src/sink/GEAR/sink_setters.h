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
#ifndef SWIFT_GEAR_STAR_FORMATION_SETTERS_H
#define SWIFT_GEAR_STAR_FORMATION_SETTERS_H

#include "sink_part.h"

/**
 * @file src/sink/GEAR/sink_setters.h
 * @brief Setters functions for GEAR star formation scheme to avoid exposing
 * implementation details to the outer world. Keep the code clean and lean.
 */


/**
 * @brief Set the birth time/scale-factor of a star particle.
 *
 * @param sp The #spart.
 * @param birth_time Birth time of the star.
 * @param birth_scale_factor Birth scale-factor of the star.
 * @param with_cosmology If we run with cosmology.
 */

__attribute__((always_inline)) INLINE void
sink_set_sink_birth_time_or_scale_factor(
    struct sink *restrict sink, const float birth_time,
    const float birth_scale_factor, const int with_cosmology) {
  if (with_cosmology) {
    sink->birth_scale_factor = birth_scale_factor;
  } else {
    sink->birth_time = birth_time;
  }
}
#endif /* SWIFT_GEAR_STAR_FORMATION_SETTERS_H */
