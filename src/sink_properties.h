/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Loic Hausammann (loic.hausammann@epfl.ch)
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
#ifndef SWIFT_SINK_PROPERTIES_H
#define SWIFT_SINK_PROPERTIES_H

/* Config parameters. */
#include <config.h>

/* Local headers. */
#include "inline.h"
#include "periodic.h"

/* Select the correct sink model */
#if defined(SINK_NONE)
#include "./sink/Default/sink_properties.h"
#elif defined(SINK_BASIC)
#include "./sink/Basic/sink_properties.h"
#elif defined(SINK_GEAR)
#include "./sink/GEAR/sink_properties.h"
#else
#error "Invalid choice of sink model"
#endif

/**
 * @brief Is the fixed-aperture gas-gas sink-formation preparation loop
 * (task subtype sink_formation_gas) active?
 *
 * The loop only makes sense when the sink model uses a single, fixed,
 * global aperture radius (GEAR's cut_off_radius). With a variable/h-based
 * cutoff there is no single radius to make the task graph respect, so the
 * loop must not be used and sink formation keeps its original path.
 *
 * This is the single gate for all formation_gas-related code: task
 * creation, dependency wiring, barrier creation, unskip activation, and the
 * various task-count/skip bookkeeping in the engine.
 *
 * @param sink_properties The #sink_props of this run, or NULL.
 */
__attribute__((always_inline)) INLINE static int
sink_formation_gas_loop_is_active(const struct sink_props *sink_properties) {
#ifdef SINKS_WITH_FIXED_CUTOFF_RADIUS
  return (sink_properties != NULL) && (sink_properties->use_fixed_r_cut != 0);
#else
  return 0;
#endif
}

/**
 * @brief Return the fixed aperture radius used by the sink-formation gas-gas
 * preparation loop.
 *
 * Only meaningful once sink_formation_gas_loop_is_active() has returned a
 * true value for the same #sink_props; returns a negative value otherwise
 * (mirrors the convention already used in runner_main.c).
 *
 * @param sink_properties The #sink_props of this run.
 */
__attribute__((always_inline)) INLINE static float
sink_formation_gas_loop_r_cut(const struct sink_props *sink_properties) {
#ifdef SINKS_WITH_FIXED_CUTOFF_RADIUS
  return sink_properties->cut_off_radius;
#else
  return -1.f;
#endif
}

/**
 * @brief Squared separation of two gas particles in the sink-formation
 * gas-gas preparation loop.
 *
 * The neighbour search and its brute-force count check both use it, so an
 * edge particle is classified the same way by both. It works in double on
 * the absolute positions, so it does not depend on the visiting cell pair.
 *
 * @param xi Position of the first particle.
 * @param xj Position of the second particle.
 * @param dim Size of the simulation box.
 * @param periodic Is the box periodic?
 * @param dx (return) Separation xi - xj of the nearest image.
 * @return The squared separation.
 */
__attribute__((always_inline)) INLINE static double sink_formation_gas_loop_r2(
    const double xi[3], const double xj[3], const double dim[3],
    const int periodic, float dx[3]) {

  double d[3] = {xi[0] - xj[0], xi[1] - xj[1], xi[2] - xj[2]};
  if (periodic) {
    d[0] = nearest(d[0], dim[0]);
    d[1] = nearest(d[1], dim[1]);
    d[2] = nearest(d[2], dim[2]);
  }
  dx[0] = (float)d[0];
  dx[1] = (float)d[1];
  dx[2] = (float)d[2];
  return d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
}

/**
 * @brief Does a squared separation from sink_formation_gas_loop_r2() lie
 * inside the aperture of the sink-formation gas-gas preparation loop?
 *
 * @param r2 The squared separation.
 * @param r_cut The aperture radius.
 */
__attribute__((always_inline)) INLINE static int
sink_formation_gas_loop_in_aperture(const double r2, const float r_cut) {
  return r2 < (double)r_cut * (double)r_cut;
}

#endif /* SWIFT_SINK_PROPERTIES_H */
