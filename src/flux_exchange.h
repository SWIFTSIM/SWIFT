/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@alumni.epfl.ch)
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
#ifndef SWIFT_FLUX_EXCHANGE_H
#define SWIFT_FLUX_EXCHANGE_H

/* Config parameters. */
#include "inline.h"

#include <config.h>

/**
 * @file flux_exchange.h
 * @brief Primitives for picking, deterministically and without shared
 * state, which of two independent visits to an active pair applies a
 * symmetric flux exactly once. Used by any module needing that guarantee
 * (see chemistry_iact.h's flux-exchange dispatcher for a worked example).
 */

/**
 * @brief Strict total order on absolute particle positions.
 *
 * Used as a tie-break so two independent evaluations of the same pair
 * agree on which particle comes first, with no shared state. Equal
 * positions mean the particles are coincident.
 *
 * Must be called with p->x (absolute, lab-frame), never a cell-relative
 * or sid-shifted local coordinate: those depend on which cell pair is
 * being processed, so the two evaluations could disagree, and a foreign
 * particle's proxy copy only matches its owner in the absolute frame.
 */
__attribute__((always_inline)) INLINE static int flux_exchange_precedes(
    const double xa[3], const double xb[3]) {
  if (xa[0] != xb[0]) return xa[0] < xb[0];
  if (xa[1] != xb[1]) return xa[1] < xb[1];
  return xa[2] < xb[2];
}

/**
 * @brief Designate the owner of a mixed-depth-band pair: the larger-h
 * particle, with flux_exchange_precedes() as an h-tie fallback.
 *
 * A pair interacts only if r2 < gamma^2 * max(ha, hb)^2, so the partner
 * always lies inside the larger-h particle's own kernel. SWIFT's neighbour
 * search therefore guarantees the larger-h particle visits this pair,
 * which a position-only tie-break cannot guarantee.
 */
__attribute__((always_inline)) INLINE static int flux_exchange_is_owner(
    const float ha, const float hb, const double xa[3], const double xb[3]) {
  if (ha != hb) return ha > hb;
  return flux_exchange_precedes(xa, xb);
}

#endif /* SWIFT_FLUX_EXCHANGE_H */
