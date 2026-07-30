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
#ifndef SWIFT_CHEMISTRY_GEAR_FVPM_DIFFUSION_DEBUG_H
#define SWIFT_CHEMISTRY_GEAR_FVPM_DIFFUSION_DEBUG_H

#include "chemistry_properties.h"
#include "error.h"
#include "hydro.h"

__attribute__((always_inline)) INLINE static void chemistry_debug_particle(
    const struct part *p, const struct xpart *xp) {

  warning("[PID%lld] chemistry_part_data:", p->id);
  for (int i = 0; i < GEAR_CHEMISTRY_ELEMENT_COUNT; i++) {
    warning("[PID%lld] metal_mass[%i]: %.3e", p->id, i,
            p->chemistry_data.metal_mass[i]);
  }
}

__attribute__((always_inline)) INLINE static void chemistry_debug_print_matrix(
    const double K[3][3]) {
  message("K = [[%e %e %e], [%e %e %e], [%e %e %e]]", K[0][0], K[0][1], K[0][2],
          K[1][0], K[1][1], K[1][2], K[2][0], K[2][1], K[2][2]);
}

#ifdef GEAR_FVPM_DIFF_DEBUG_PAIR_VISIT_COUNT

/* Per-pair visit counters for the flux-exchange dispatcher (chemistry_iact.h).
   Aggregate acted/skipped alone cannot prove no pair was dropped (a
   both_updatable_here visit covers a pair without touching either counter),
   so the per-pair log below is the actual proof; these are a cheap sanity
   readout. Defined in chemistry_debug.c. */
extern long long chemistry_fvpm_visit_acted;
extern long long chemistry_fvpm_visit_skipped;

/* Per-pair log keyed by particle address, used to verify every pair is
   covered exactly once. Not thread-safe: diagnostic use only,
   --threads=1. Defined in chemistry_debug.c. */
extern void chemistry_fvpm_visit_log_record(const struct part *pi,
                                            const struct part *pj,
                                            enum chemistry_fvpm_visit_arm arm);
extern void chemistry_fvpm_visit_log_analyze(void);
extern void chemistry_fvpm_visit_log_reset(void);

/* atexit() callback: reports the acted/skipped totals and the per-pair
   log's dropped/double_applied verdict. Defined in chemistry_debug.c. */
extern void chemistry_fvpm_print_visit_counts(void);

#endif /* GEAR_FVPM_DIFF_DEBUG_PAIR_VISIT_COUNT */

#endif /* SWIFT_CHEMISTRY_GEAR_DEBUG_H */
