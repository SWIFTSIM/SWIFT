/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2015 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2016 Tom Theuns (tom.theuns@durham.ac.uk)
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

#include <float.h>
#include "feedback.h"

/**
 * @brief Computes the feedback time-step of a given particle due to
 * a source term
 *
 * This function only branches towards the potential chosen by the user.
 *
 * @param feedback The properties of the source terms.
 * @param phys_const The physical constants in internal units.
 * @param g Pointer to the particle data.
 */
__attribute__((always_inline)) INLINE static float sourceterms_compute_timestep(
    const struct sourceterms* feedback,
    const struct phys_const* const phys_const, const struct part* const p) {
  float dt = FLT_MAX;
#ifdef SN_FEEDBACK
  dt = fmin(dt, sn_feedback_timestep(feedback, phys_const, p));
#endif
}
__attribute__((always_inline)) INLINE static void feedback(
    const struct sourceterms* feedback,
    const struct phys_const* const phys_const, struct part* p)
#ifdef SN_FEEDBACK
    sn_feedback(feedback, phys_const, p);
#endif
}
