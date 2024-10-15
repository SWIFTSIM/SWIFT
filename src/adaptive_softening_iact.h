/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2023 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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

/* Config parameters. */
#ifndef SWIFT_ADAPTIVE_SOFTENING_IACT_H
#define SWIFT_ADAPTIVE_SOFTENING_IACT_H

/* Config parameters. */
#include <config.h>

/* Local headers. */
#include "adaptive_softening_struct.h"
#include "inline.h"
#include "kernel_hydro.h"

#ifdef ADAPTIVE_SOFTENING

/**
 * @brief Computes the contribution to the softening length change term.
 *
 * Computes equation 26 of Price & Monaghan 2007, MNRAS, 374, 4.
 * without the prefactor. The prefactor gets added in end_density().
 *
 * @param pi The #part for which we compute terms.
 * @param ui The ratio of the inter-particle distance to the smoothing length.
 * @param hi_inv The inverse the particle's smoothing length.
 * @param mj The mass of the other particle.
 */
__attribute__((always_inline)) INLINE static void
adaptive_softening_add_correction_term(struct part *pi, const float ui,
                                       const float hi_inv, const float mj) {

  pi->adaptive_softening_data.zeta += mj * potential_dh(ui, hi_inv);
}

/**
 * @brief Computes the norm of the accleration due to the change in softening.
 *
 * Computes equation 27 of Price & Monaghan 2007, MNRAS, 374, 4.
 * without the mass and the distance vector. These are multiplied in
 * by the parent function.
 *
 * @param pi The first #part.
 * @param pj The second #part.
 * @param wi_dr The norm of the kernel gradient for i.
 * @param wj_dr The norm of the kernel gradient for j.
 * @param f_ij The smoothling-length correction term for i.
 * @param f_ji The smoothling-length correction term for j.
 * @param r_inv the inverse of the distance linking the particles.
 */
__attribute__((always_inline)) INLINE static float
adaptive_softening_get_acc_term(const struct part *restrict pi,
                                const struct part *restrict pj,
                                const float wi_dr, const float wj_dr,
                                const float f_ij, const float f_ji,
                                const float r_inv) {
  /* Recover some data */
  const float zetai = pi->adaptive_softening_data.zeta;
  const float zetaj = pj->adaptive_softening_data.zeta;

  /* Adaptive softening acceleration term
   * Price & Monaghan 2007, eq. 27 (second term)
   * Note that G/2 is included in the zeta terms.
   * Note also that f_ij is 1/Omega_i in Price's notation. */
  const float adap_acc_term =
      0.5f * (zetai * f_ij * wi_dr + zetaj * f_ji * wj_dr) * r_inv;

  return adap_acc_term;
}

#else

/**
 * @brief Computes the contribution to the softening length change term.
 *
 * No adaptive softening --> Nothing to do.
 *
 * @param pi The #part for which we compute terms.
 * @param ui The ratio of the inter-particle distance to the smoothing length.
 * @param hi_inv The inverse the particle's smoothing length.
 * @param mj The mass of the other particle.
 */
__attribute__((always_inline)) INLINE static void
adaptive_softening_add_correction_term(struct part *pi, const float ui,
                                       const float hi_inv, const float mj) {}

/**
 * @brief Computes the norm of the accleration due to the change in softening.
 *
 * No adaptive softening --> Nothing to do --> Return 0.
 *
 * @param pi The first #part.
 * @param pj The second #part.
 * @param wi_dr The norm of the kernel gradient for i.
 * @param wj_dr The norm of the kernel gradient for j.
 * @param f_ij The smoothling-length correction term for i.
 * @param f_ji The smoothling-length correction term for j.
 * @param r_inv the inverse of the distance linking the particles.
 */
__attribute__((always_inline)) INLINE static float
adaptive_softening_get_acc_term(const struct part *restrict pi,
                                const struct part *restrict pj,
                                const float wi_dr, const float wj_dr,
                                const float f_ij, const float f_ji,
                                const float r_inv) {
  return 0.f;
}

#endif /* ADAPTIVE_SOFTENING */

#endif /* SWIFT_ADAPTIVE_SOFTENING_IACT_H */
