/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2022 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_TRACERS_EAGLE_DEBUG_H
#define SWIFT_TRACERS_EAGLE_DEBUG_H

__attribute__((always_inline)) INLINE static void tracers_debug_particle(
    const struct part* p, const struct xpart* xp) {

  if (xp != NULL) {
    warning("[PID%lld] tracers_xpart_data:", p->id);
    warning(
        "[PID%lld] maximum_temperature = %.3e, "
        "maximum_temperature_scale_factor/time = %.3e, "
        "last_AGN_injection_scale_factor/time = %.3e, "
        "density_before_last_AGN_feedback_event = %.3e, "
        "entropy_before_last_AGN_feedback_event = %.3e, "
        "density_at_last_AGN_feedback_event = %.3e, "
        "entropy_at_last_AGN_feedback_event = %.3e, AGN_feedback_energy = "
        "%.3e, "
        "hit_by_SNII_feedback = %d, hit_by_AGN_feedback = %d",
        p->id, xp->tracers_data.maximum_temperature,
        xp->tracers_data.maximum_temperature_scale_factor,
        xp->tracers_data.last_AGN_injection_scale_factor,
        xp->tracers_data.density_before_last_AGN_feedback_event,
        xp->tracers_data.entropy_before_last_AGN_feedback_event,
        xp->tracers_data.density_at_last_AGN_feedback_event,
        xp->tracers_data.entropy_at_last_AGN_feedback_event,
        xp->tracers_data.AGN_feedback_energy,
        xp->tracers_data.hit_by_SNII_feedback,
        xp->tracers_data.hit_by_AGN_feedback);
  }
}

#endif /* SWIFT_TRACERS_EAGLE_DEBUG_H */
