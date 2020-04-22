/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#ifndef SWIFT_MULTIPOLE_ACCEPT_H
#define SWIFT_MULTIPOLE_ACCEPT_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "multipole_struct.h"

/**
 * @brief Checks whether a cell-cell interaction can be appromixated by a M-M
 * interaction using the distance and cell radius.
 *
 * We use the multipole acceptance criterion of Dehnen, 2002, JCoPh, Volume 179,
 * Issue 1, pp.27-42, equation 10.
 *
 * We also additionally check that the distance between the multipoles
 * is larger than the softening lengths (here the distance at which
 * the gravity becomes Newtonian again, not the Plummer-equivalent quantity).
 *
 * @param r_crit_a The size of the multipole A.
 * @param r_crit_b The size of the multipole B.
 * @param theta_crit2 The square of the critical opening angle.
 * @param r2 Square of the distance (periodically wrapped) between the
 * multipoles.
 * @param epsilon_a The maximal softening length of any particle in A.
 * @param epsilon_b The maximal softening length of any particle in B.
 */
__attribute__((always_inline, const)) INLINE static int gravity_M2L_accept(
    const double r_crit_a, const double r_crit_b, const double theta_crit2,
    const double r2, const double epsilon_a, const double epsilon_b) {

  const double size = r_crit_a + r_crit_b;
  const double size2 = size * size;
  const double epsilon_a2 = epsilon_a * epsilon_a;
  const double epsilon_b2 = epsilon_b * epsilon_b;

  // MATTHIEU: Make this mass-dependent ?

  /* Multipole acceptance criterion (Dehnen 2002, eq.10) */
  return (r2 * theta_crit2 > size2) && (r2 > epsilon_a2) && (r2 > epsilon_b2);
}

/**
 * @brief Checks whether a particle-cell interaction can be appromixated by a
 * M2P interaction using the distance and cell radius.
 *
 * We use the multipole acceptance criterion of Dehnen, 2002, JCoPh, Volume 179,
 * Issue 1, pp.27-42, equation 10.
 *
 * We also additionally check that the distance between the particle and the
 * multipole is larger than the softening length (here the distance at which
 * the gravity becomes Newtonian again, not the Plummer-equivalent quantity).
 *
 * @param r_max2 The square of the size of the multipole.
 * @param theta_crit2 The square of the critical opening angle.
 * @param r2 Square of the distance (periodically wrapped) between the
 * particle and the multipole.
 * @param epsilon The softening length of the particle.
 */
__attribute__((always_inline, const)) INLINE static int gravity_M2P_accept(
    const float r_max2, const float theta_crit2, const float r2,
    const float epsilon) {

  // MATTHIEU: Make this mass-dependent ?

  /* Multipole acceptance criterion (Dehnen 2002, eq.10) */
  return (r2 * theta_crit2 > r_max2) && (r2 > epsilon * epsilon);
}

#endif /* SWIFT_MULTIPOLE_ACCEPT_H */
