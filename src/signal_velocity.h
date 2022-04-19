/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2022 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_SIGNAL_VELOCITY_H
#define SWIFT_SIGNAL_VELOCITY_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "part.h"

#ifndef NONE_MHD

/**
 * @brief Compute the signal velocity between two gas particles.
 *
 * MHD case.
 * This is eq. (131) of Price D., JCoPh, 2012, Vol. 231, Issue 3.
 *
 * @brief pi The first #part.
 * @brief pj The second #part.
 * @brief mu_ij The velocity on the axis linking the particles, or zero if the
 * particles are moving away from each other,
 * @brief beta The non-linear viscosity constant.
 */
__attribute__((always_inline)) INLINE static float signal_velocity(
    const struct part *restrict pi, const struct part *restrict pj,
    const float mu_ij, const float beta) {

  // TODO: Implement !

  return -1.f;
}

#else

/**
 * @brief Compute the signal velocity between two gas particles.
 *
 * Non-MHD case.
 * This is eq. (103) of Price D., JCoPh, 2012, Vol. 231, Issue 3.
 *
 * @brief pi The first #part.
 * @brief pj The second #part.
 * @brief mu_ij The velocity on the axis linking the particles, or zero if the
 * particles are moving away from each other,
 * @brief beta The non-linear viscosity constant.
 */
__attribute__((always_inline)) INLINE static float signal_velocity(
    const struct part *restrict pi, const struct part *restrict pj,
    const float mu_ij, const float beta) {

  const float ci = pi->force.soundspeed;
  const float cj = pj->force.soundspeed;

  return ci + cj - beta * mu_ij;
}

#endif

#endif /* SWIFT_SIGNAL_VELOCITY_H */
