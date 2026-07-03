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

#endif /* SWIFT_SINK_PROPERTIES_H */
