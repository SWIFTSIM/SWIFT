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
#include <config.h>

#ifndef NONE_MHD

/* Include MHD definition of signal velocity */
#include "mhd.h"

/**
 * @brief Compute the signal velocity between two gas particles,
 * MHD case.
 *
 * Warning ONLY to be called just after preparation of the force loop.
 *
 * @param dx Comoving vector separating both particles (pi - pj).
 * @brief pi The first #part.
 * @brief pj The second #part.
 * @brief mu_ij The velocity on the axis linking the particles, or zero if the
 * particles are moving away from each other,
 * @brief beta The non-linear viscosity constant.
 */
__attribute__((always_inline)) INLINE static float signal_velocity(
    const float dx[3], const struct part *restrict pi,
    const struct part *restrict pj, const float mu_ij, const float beta) {

  return mhd_signal_velocity(dx, pi, pj, mu_ij, beta);
}

#else

/* Include hydro definition of signal velocity */
#include "hydro.h"

/**
 * @brief Compute the signal velocity between two gas particles,
 * Non-MHD case.
 *
 * Warning: Can ONLY to be called just after preparation of the force loop.
 *
 * @param dx Comoving vector separating both particles (pi - pj).
 * @brief pi The first #part.
 * @brief pj The second #part.
 * @brief mu_ij The velocity on the axis linking the particles, or zero if the
 * particles are moving away from each other,
 * @brief beta The non-linear viscosity constant.
 */
__attribute__((always_inline)) INLINE static float signal_velocity(
    const float dx[3], const struct part *restrict pi,
    const struct part *restrict pj, const float mu_ij, const float beta) {

  return hydro_signal_velocity(dx, pi, pj, mu_ij, beta);
}

#endif

#endif /* SWIFT_SIGNAL_VELOCITY_H */
