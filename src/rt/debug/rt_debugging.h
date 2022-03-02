/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2021 Mladen Ivkovic (mladen.ivkovic@hotmail.com)
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
#ifndef SWIFT_RT_DEBUGGING_DEBUG_H
#define SWIFT_RT_DEBUGGING_DEBUG_H

#include "rt_properties.h"

/**
 * @file src/rt/debug/rt_debugging.h
 * @brief Main header file for the debug radiative transfer scheme
 * extra debugging functions.
 */

/**
 * @brief Debugging checks loop over all star particles after each time step
 */
static void rt_debugging_end_of_step_stars_mapper(void *restrict map_data,
                                                  int scount,
                                                  void *restrict extra_data) {

  struct spart *restrict sparts = (struct spart *)map_data;
  const struct engine *restrict e = (struct engine *)extra_data;

  unsigned long long emission_sum_this_step = 0ULL;
  unsigned long long emission_sum_tot = 0ULL;

  for (int k = 0; k < scount; k++) {

    struct spart *restrict sp = &sparts[k];
    emission_sum_this_step += sp->rt_data.debug_iact_hydro_inject;
    emission_sum_tot += sp->rt_data.debug_radiation_emitted_tot;
    /* Reset all values here in case stars won't be active next step */
    sp->rt_data.debug_iact_hydro_inject = 0;
    sp->rt_data.debug_iact_hydro_inject_prep = 0;
  }

  atomic_add(&e->rt_props->debug_radiation_emitted_this_step,
             emission_sum_this_step);
  atomic_add(&e->rt_props->debug_radiation_emitted_tot, emission_sum_tot);
}

/**
 * @brief Debugging checks loop over all hydro particles after each time step
 */
static void rt_debugging_end_of_step_hydro_mapper(void *restrict map_data,
                                                  int count,
                                                  void *restrict extra_data) {

  struct part *restrict parts = (struct part *)map_data;
  const struct engine *restrict e = (struct engine *)extra_data;

  unsigned long long absorption_sum_this_step = 0ULL;
  unsigned long long absorption_sum_tot = 0ULL;

  for (int k = 0; k < count; k++) {

    struct part *restrict p = &parts[k];
    absorption_sum_this_step += p->rt_data.debug_iact_stars_inject;
    absorption_sum_tot += p->rt_data.debug_radiation_absorbed_tot;
    /* Reset all values here in case particles won't be active next step */
    p->rt_data.debug_iact_stars_inject = 0;
  }

  atomic_add(&e->rt_props->debug_radiation_absorbed_this_step,
             absorption_sum_this_step);
  atomic_add(&e->rt_props->debug_radiation_absorbed_tot, absorption_sum_tot);
}

/**
 * @brief At the end of each time step, loop over both hydro and star
 * particles and do whatever checks for this particular time step you
 * want done.
 *
 * @param e The #engine.
 * @param verbose Are we talkative?
 */
__attribute__((always_inline)) INLINE static void
rt_debugging_checks_end_of_step(struct engine *e, int verbose) {

  struct space *s = e->s;
  if (!(e->policy & engine_policy_rt)) return;
#ifdef WITH_MPI
  /* Since we aren't sending data back, none of these checks will
   * pass a run over MPI. */
  return;
#endif

  const ticks tic = getticks();

  /* reset values before the particle loops.
   * reset total counts as well. We track the totals since the beginning
   * of time in particles individually. */
  e->rt_props->debug_radiation_emitted_this_step = 0ULL;
  e->rt_props->debug_radiation_absorbed_this_step = 0ULL;
  e->rt_props->debug_radiation_emitted_tot = 0ULL;
  e->rt_props->debug_radiation_absorbed_tot = 0ULL;

  /* hydro particle loop */
  if (s->nr_parts > 0)
    threadpool_map(&e->threadpool, rt_debugging_end_of_step_hydro_mapper,
                   s->parts, s->nr_parts, sizeof(struct part),
                   threadpool_auto_chunk_size, /*extra_data=*/e);

  /* star particle loop */
  if (s->nr_sparts > 0)
    threadpool_map(&e->threadpool, rt_debugging_end_of_step_stars_mapper,
                   s->sparts, s->nr_sparts, sizeof(struct spart),
                   threadpool_auto_chunk_size, /*extra_data=*/e);

  /* Have we accidentally invented or deleted some radiation somewhere? */
  if ((e->rt_props->debug_radiation_emitted_this_step !=
       e->rt_props->debug_radiation_absorbed_this_step) ||
      (e->rt_props->debug_radiation_emitted_tot !=
       e->rt_props->debug_radiation_absorbed_tot))
    error(
        "Emitted and absorbed radiation vary.\n"
        "  This step: star emission %12lld; gas absorption %12lld\n"
        "Since start: star emission %12lld; gas absorption %12lld",
        e->rt_props->debug_radiation_emitted_this_step,
        e->rt_props->debug_radiation_absorbed_this_step,
        e->rt_props->debug_radiation_emitted_tot,
        e->rt_props->debug_radiation_absorbed_tot);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

#endif /* SWIFT_RT_DEBUGGING_DEBUG_H */
