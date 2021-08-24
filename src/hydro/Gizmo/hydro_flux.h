/*******************************************************************************
 * This file is part of SWIFT.
 * Coypright (c) 2019 Bert Vandenbroucke (bert.vandenbroucke@gmail.com)
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
#ifndef SWIFT_GIZMO_HYDRO_FLUX_H
#define SWIFT_GIZMO_HYDRO_FLUX_H

#if defined(GIZMO_MFV_SPH)
#include "MFV/hydro_flux.h"
#elif defined(GIZMO_MFM_SPH)
#include "MFM/hydro_flux.h"
#endif

__attribute__((always_inline)) INLINE static void hydro_part_reset_flux_count(
    struct part* restrict p) {
  p->flux.flux_count = 0;
}

__attribute__((always_inline)) INLINE static void
hydro_part_update_flux_count_left(struct part* restrict pi,
                                  const struct part* restrict pj) {
  if (pi->id < pj->id) {
    pi->flux.flux_count += 1;
  } else if (pj->id < pi->id) {
    pi->flux.flux_count -= 1;
  } else {
    error("Particle IDs are the same!");
  }
}

__attribute__((always_inline)) INLINE static void
hydro_part_update_flux_count_right(const struct part* restrict pi,
                                   struct part* restrict pj) {
  if (pi->id < pj->id) {
    pj->flux.flux_count -= 1;
  } else if (pj->id < pi->id) {
    pj->flux.flux_count += 1;
  } else {
    error("Particle IDs are the same!");
  }
}

#endif /* SWIFT_GIZMO_HYDRO_FLUX_H */
