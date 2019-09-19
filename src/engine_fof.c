/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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
#include "../config.h"

/* This object's header. */
#include "engine.h"

/**
 * @brief Activate all the #gpart communications in preparation
 * fof a call to FOF.
 *
 * @param e The #engine to act on.
 */
void engine_activate_gpart_comms(struct engine *e) {

#ifdef WITH_MPI

  const ticks tic = getticks();

  struct scheduler *s = &e->sched;
  const int nr_tasks = s->nr_tasks;
  struct task *tasks = s->tasks;

  for (int k = 0; k < nr_tasks; ++k) {

    struct task *t = &tasks[k];

    if ((t->type == task_type_send) && (t->subtype == task_subtype_gpart)) {
      scheduler_activate(s, t);
    } else if ((t->type == task_type_recv) &&
               (t->subtype == task_subtype_gpart)) {
      scheduler_activate(s, t);
    } else {
      t->skip = 1;
    }
  }

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());

#else
  error("Calling an MPI function in non-MPI mode.");
#endif
}

/**
 * @brief Activate all the FOF tasks.
 *
 * Marks all the other task types to be skipped.
 *
 * @param e The #engine to act on.
 */
void engine_activate_fof_tasks(struct engine *e) {

  const ticks tic = getticks();

  struct scheduler *s = &e->sched;
  const int nr_tasks = s->nr_tasks;
  struct task *tasks = s->tasks;

  for (int k = 0; k < nr_tasks; k++) {

    struct task *t = &tasks[k];

    if (t->type == task_type_fof_self || t->type == task_type_fof_pair)
      scheduler_activate(s, t);
    else
      t->skip = 1;
  }

  if (e->verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

/**
 * @brief Run a FOF search.
 *
 * @param e the engine
 * @param dump_results Are we writing group catalogues to output files?
 * @param seed_black_holes Are we seeding black holes?
 */
void engine_fof(struct engine *e, const int dump_results,
                const int seed_black_holes) {

#ifdef WITH_FOF

  ticks tic = getticks();

  /* Compute number of DM particles */
  const long long total_nr_baryons =
      e->total_nr_parts + e->total_nr_sparts + e->total_nr_bparts;
  const long long total_nr_dmparts =
      e->total_nr_gparts - e->total_nr_DM_background_gparts - total_nr_baryons;

  /* Initialise FOF parameters and allocate FOF arrays. */
  fof_allocate(e->s, total_nr_dmparts, e->fof_properties);

  /* Make FOF tasks */
  engine_make_fof_tasks(e);

  /* and activate them. */
  engine_activate_fof_tasks(e);

  /* Perform local FOF tasks. */
  engine_launch(e);

  /* Perform FOF search over foreign particles and
   * find groups which require black hole seeding.  */
  fof_search_tree(e->fof_properties, e->black_holes_properties,
                  e->physical_constants, e->cosmology, e->s, dump_results,
                  seed_black_holes);

  /* Reset flag. */
  e->run_fof = 0;

  /* Flag that a FOF has taken place */
  e->step_props |= engine_step_prop_fof;

  /* ... and find the next FOF time */
  if (seed_black_holes) engine_compute_next_fof_time(e);

  if (engine_rank == 0)
    message("Complete FOF search took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());
#else
  error("SWIFT was not compiled with FOF enabled!");
#endif
}
