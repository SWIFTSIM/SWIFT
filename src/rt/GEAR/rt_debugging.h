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
#ifndef SWIFT_RT_DEBUGGING_GEAR_H
#define SWIFT_RT_DEBUGGING_GEAR_H

#ifdef SWIFT_RT_DEBUG_CHECKS

#include "active.h"
#include "rt_properties.h"
#include "timeline.h"

/**
 * @file src/rt/GEAR/rt_debugging.h
 * @brief Main header file for the GEAR radiative transfer scheme
 * extra debugging functions.
 */

/**
 * @brief Check whether RT time step size is valid compared to the hydro one.
 *
 * @param p particle to work on
 * @param dti_rt current RT integer time step
 * @param dti_hydro current hydro integer time step
 * @param max_nr_rt_subcycles max number of RT sub-cycles.
 * @param time_base minimal time step in this run
 */
__attribute__((always_inline)) INLINE static void rt_debugging_check_timestep(
    const struct part *restrict p, integertime_t *dti_rt,
    const integertime_t *dti_hydro, int max_nr_rt_subcycles, double time_base) {

  if (*dti_hydro < *dti_rt)
    error("part %lld has hydro time (%lld/%g) < RT time step (%lld/%g", p->id,
          *dti_hydro, *dti_hydro * time_base, *dti_rt, *dti_rt * time_base);
}

/**
 * @brief This resets particle carried quantities after each subcycling
 * step such that the internal checks are still consistent.
 * @param p the particle to work on
 */
__attribute__((always_inline)) INLINE static void rt_debugging_count_subcycle(
    struct part *restrict p) {
  p->rt_data.debug_nsubcycles += 1;
}

/**
 * @brief Check that the particle completed the correct number of subcycles.
 * This is checked in every rt_reset_part, before the subcycling count is reset.
 * @param p the particle to work on
 * @param rt_props RT properties struct
 */
__attribute__((always_inline)) INLINE static void
rt_debugging_check_nr_subcycles(struct part *restrict p,
                                const struct rt_props *rt_props) {

  /* TODO: this check may fail when running with limiter/sync. */

  /* NOTE: we need to do this check somewhere in the hydro tasks.
   * (1) it needs to be done when all tasks are active and before the
   * particle hydro time step is changed.
   * (2) If you do it during RT tasks, it won't properly check how
   * many sub-cycles you did during a single hydro task.
   * (3) You can't do it during the timestep task, since between
   * the hydro and the timestep we already do an RT step. */

  /* skip initialization */
  if (p->time_bin == 0) return;
  if (p->rt_time_data.time_bin == 0)
    error("Got part %lld with RT time bin 0", p->id);

  timebin_t bindiff = p->time_bin - p->rt_time_data.time_bin;

  if (rt_props->debug_max_nr_subcycles <= 1) {
    /* Running without subcycling. */
    if (bindiff != 0)
      error("Running without subcycling but got bindiff=%d for part=%lld",
            bindiff, p->id);
    if (p->rt_data.debug_nsubcycles != 1)
      error("Running without subcycling but got part=%lld subcycle count=%d",
            p->id, p->rt_data.debug_nsubcycles);
    return;
  }
}

/**
 * @brief This resets particle carried quantities after each subcycling
 * step such that the internal checks are still consistent.
 * @param p the particle to work on
 */
__attribute__((always_inline)) INLINE static void
rt_debugging_reset_each_subcycle(struct part *restrict p) {

  p->rt_data.debug_calls_iact_gradient_interaction = 0;
  p->rt_data.debug_calls_iact_transport_interaction = 0;

  p->rt_data.debug_injection_done = 0;
  p->rt_data.debug_gradients_done = 0;
  p->rt_data.debug_transport_done = 0;
  p->rt_data.debug_thermochem_done = 0;
}

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

#ifndef WITH_MPI
    /* These checks will fail with MPI since we're not
     * sending data back */
    for (int g = 0; g < RT_NGROUPS; g++) {
      /* also check now that we actually injected the correct
       * amount of energy
       * sp->rt_data.emission_this_step: energy we should distribute
       *                                 this step
       * sp->rt_data.debug_injected_energy: energy we actually did
       *                                    distribute this step    */
      if (sp->rt_data.debug_injected_energy[g] != 0.f) {
        float diff = 1.f - sp->rt_data.emission_this_step[g] /
                               sp->rt_data.debug_injected_energy[g];

        if (fabsf(diff) > 1e-4) {
          /* Dividing the total into several parts and summing them up again
           * while hoping to obtain the same results may lead to diappointment
           * due to roundoff errors. Check that the sum of the individual
           * weights and the ones we collected for the injection are close
           * enough. */
          float psi_sum_now = 0.f;
          for (int i = 0; i < 8; i++)
            psi_sum_now += sp->rt_data.octant_weights[i];
          float diff_weights = 1.f - sp->rt_data.debug_psi_sum / psi_sum_now;
          if (fabsf(diff_weights) > 1e-4)
            message(
                "Incorrect injection ID %lld: "
                "group %d expected %.3g got %.3g diff %.3g diff_weights %.3g",
                sp->id, g, sp->rt_data.emission_this_step[g],
                sp->rt_data.debug_injected_energy[g], diff, diff_weights);
        }
      }
    }
#endif

    for (int g = 0; g < RT_NGROUPS; g++) {
      sp->rt_data.debug_injected_energy[g] = 0.f;
    }
    for (int g = 0; g < RT_NGROUPS; g++) {
      sp->rt_data.emission_this_step[g] = 0.f;
    }
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
  float energy_sum[RT_NGROUPS];
  for (int g = 0; g < RT_NGROUPS; g++) energy_sum[g] = 0.f;

  for (int k = 0; k < count; k++) {

    struct part *restrict p = &parts[k];
    absorption_sum_this_step += p->rt_data.debug_iact_stars_inject;
    absorption_sum_tot += p->rt_data.debug_radiation_absorbed_tot;

    /* Reset all values here in case particles won't be active next step */
    p->rt_data.debug_iact_stars_inject = 0;

    /* Sum up total energies for budget */
    for (int g = 0; g < RT_NGROUPS; g++) {
      energy_sum[g] +=
          p->rt_data.radiation[g].energy_density * p->geometry.volume;
    }
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

#ifdef WITH_MPI
  /* Since we aren't sending data back, none of these checks will
   * pass a run over MPI. Make sure you run the threadpool functions
   * first though so certain variables can get reset properly. */
  return;
#endif

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

/**
 * @brief Perform a series of consistency and sanity checks.
 *
 * @param p particle to check
 * @param loc location where this is called from. This determines which checks
 * will be done:
 *
 * 0: during kicks/after drifts.
 * 1: during rt_ghost1/finalise_injection / after kicks.
 * 2: during gradients / after injection.
 * 3: during transport / after gradients.
 * 4: during thermochem / after transport.
 * 5: after thermochem.
 *
 * @param function_name: Function name (or message) you want printed on error.
 */
__attribute__((always_inline)) INLINE static void rt_debug_sequence_check(
    struct part *restrict p, int loc, const char *function_name) {

  /* Note: Checking whether a particle has been drifted at all is not
   * compatible with subcycling. There is no reliable point where to
   * reset the counters and have sensible results. */

  if (loc > 0) {
    /* Are kicks done? */
    if (p->rt_data.debug_nsubcycles == 0) {
      if (p->rt_data.debug_kicked != 1)
        error(
            "called %s on particle %lld with wrong kick count=%d (expected "
            "1) cycle=%d",
            function_name, p->id, p->rt_data.debug_kicked,
            p->rt_data.debug_nsubcycles);
    } else if (p->rt_data.debug_nsubcycles > 0) {
      /* This covers case 1, 2, 3, 4, 5 */
      if (p->rt_data.debug_kicked != 2)
        error(
            "called %s on particle %lld with wrong kick count=%d (expected 2) "
            "cycle=%d",
            function_name, p->id, p->rt_data.debug_kicked,
            p->rt_data.debug_nsubcycles);
    } else {
      error("Got negative subcycle???");
    }
  }

  if (loc > 1) {
    /* is injection done? */
    if (p->rt_data.debug_injection_done != 1)
      error("called %s on part %lld when finalise injection count is %d ID",
            function_name, p->id, p->rt_data.debug_injection_done);
  }

  if (loc > 2) {
    /* are gradients done? */
    if (p->rt_data.debug_gradients_done != 1)
      error("called %s on part %lld when gradients_done count is %d",
            function_name, p->id, p->rt_data.debug_gradients_done);
  }

  if (loc > 3) {
    /* is transport done? */
    if (p->rt_data.debug_transport_done != 1)
      error("called %s on part %lld when transport_done != 1: %d",
            function_name, p->id, p->rt_data.debug_transport_done);
  }

  if (loc > 4) {
    /* is thermochemistry done? */
    if (p->rt_data.debug_thermochem_done != 1)
      error("called %s on part %lld with thermochem_done count=%d",
            function_name, p->id, p->rt_data.debug_thermochem_done);
  }
}

#endif /* SWIFT_RT_DEBUG_CHECKS */
#endif /* SWIFT_RT_DEBUGGING_GEAR_H */
