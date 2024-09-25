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
#ifndef SWIFT_FEEDBACK_GEAR_MECHANICAL_DEBUG_H
#define SWIFT_FEEDBACK_GEAR_MECHANICAL_DEBUG_H

__attribute__((always_inline)) INLINE static void feedback_debug_particle(
    const struct part* p, const struct xpart* xp) {

  if (xp != NULL) {
    warning("[PID%lld] feedback_xpart_data:", p->id);
    warning(
        "[PID%lld] delta_mass = %.3e, delta_u = %.3e, delta_p = [%.3e, %.3e, "
        "%.3e]",
        p->id, xp->feedback_data.delta_mass, xp->feedback_data.delta_u,
        xp->feedback_data.delta_p[0], xp->feedback_data.delta_p[1],
        xp->feedback_data.delta_p[2]);
  }
}

#endif /* SWIFT_FEEDBACK_GEAR_MECHANICAL_DEBUG_H */
