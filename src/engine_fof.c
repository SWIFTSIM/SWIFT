/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include <config.h>

/* This object's header. */
#include "engine.h"

/* Local headers. */
#include "fof.h"

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

    if ((t->type == task_type_send) && (t->subtype == task_subtype_fof)) {
      scheduler_activate(s, t);
    } else if ((t->type == task_type_pack) &&
               (t->subtype == task_subtype_fof)) {
      scheduler_activate(s, t);
    } else if ((t->type == task_type_recv) &&
               (t->subtype == task_subtype_fof)) {
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
 * @brief Activate all the FOF linking tasks.
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
 * @brief Activate all the FOF attaching tasks.
 *
 * Marks all the other task types to be skipped.
 *
 * @param e The #engine to act on.
 */
void engine_activate_fof_attach_tasks(struct engine *e) {

  const ticks tic = getticks();

  struct scheduler *s = &e->sched;
  const int nr_tasks = s->nr_tasks;
  struct task *tasks = s->tasks;

  for (int k = 0; k < nr_tasks; k++) {

    struct task *t = &tasks[k];

    if (t->type == task_type_fof_attach_self ||
        t->type == task_type_fof_attach_pair)
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
 * @param dump_debug_results Are we writing a txt-file debug catalogue
 * (including BH seed info)?
 * @param seed_black_holes Are we seeding black holes?
 * @param foreign_buffers_allocated Are the foreign buffers currently
 * allocated?
 */
void engine_fof(struct engine *e, const int dump_results,
                const int dump_debug_results, const int seed_black_holes,
                const int foreign_buffers_allocated) {

#ifdef WITH_FOF

  const ticks tic = getticks();

  /* TEMPORARY: FOF scaling timers. */
  fof_timers_reset();
  fof_counters[fof_counter_nr_gparts].value = e->s->nr_gparts;
  fof_counters[fof_counter_nr_threads].value = e->nr_threads;
  ticks tic_phase = getticks();

  /* Start by cleaning up the foreign buffers */
  if (foreign_buffers_allocated) {
#ifdef WITH_MPI
    space_free_foreign_parts(e->s, /*clear pointers=*/1);
#endif
  }

  fof_timers[fof_timer_free_start] = getticks() - tic_phase;
  tic_phase = getticks();

  /* Free all the memory related to gas and stars sort arrays. */
  space_free_sort_indices(e->s);

  /* Initialise FOF parameters and allocate FOF arrays. */
  fof_allocate(e->s, e->fof_properties);

  fof_timers[fof_timer_allocate] = getticks() - tic_phase;
  tic_phase = getticks();

  /* Make FOF tasks */
  engine_make_fof_tasks(e);

  fof_timers[fof_timer_make_tasks] = getticks() - tic_phase;
  tic_phase = getticks();

  /* and activate them. */
  engine_activate_fof_tasks(e);

  fof_timers[fof_timer_activate_tasks] = getticks() - tic_phase;

  /* Print the number of active tasks ? */
  if (e->verbose) engine_print_task_counts(e);

  /* Perform local FOF tasks for linkable particles. */
  tic_phase = getticks();
  engine_launch(e, "fof");
  fof_timers[fof_timer_local_link] = getticks() - tic_phase;

  /* Compute group sizes (only of local fragments with MPI) */
  tic_phase = getticks();
  fof_compute_local_sizes(e->fof_properties, e->s);
  fof_timers[fof_timer_local_sizes] = getticks() - tic_phase;

#ifdef WITH_MPI

  /* Allocate buffers to receive the gpart fof information */
  tic_phase = getticks();
  engine_allocate_foreign_particles(e, /*fof=*/1);
  fof_timers[fof_timer_alloc_foreign] = getticks() - tic_phase;
  fof_counters[fof_counter_gparts_foreign_recv].value =
      e->s->nr_gparts_foreign;

  /* Compute the local<->foreign group links (nothing to do without MPI)*/
  fof_search_foreign_cells(e->fof_properties, e->s);

  /* Link the foreign fragments and finalise global group list (nothing to do
   * without MPI) */
  fof_link_foreign_fragments(e->fof_properties, e->s);
#endif

  /* Count the total number of (linkable) groups,
   * assign unique group IDs and set the IDs in the particles */
  fof_assign_group_ids(e->fof_properties, e->s);

#ifdef WITH_MPI

  /* Activate the tasks exchanging all the required gparts */
  tic_phase = getticks();
  engine_activate_gpart_comms(e);

  /* Perform send and receive tasks. */
  engine_launch(e, "fof comms");
  fof_timers[fof_timer_attach_comms] = getticks() - tic_phase;

#endif

  /* Compute the attachable->linkable links */
  tic_phase = getticks();
  fof_link_attachable_particles(e->fof_properties, e->s);
  fof_timers[fof_timer_attach] = getticks() - tic_phase;

#ifdef WITH_MPI

  /* Free the foreign particles */
  tic_phase = getticks();
  space_free_foreign_parts(e->s, /*clear pointers=*/1);
  fof_timers[fof_timer_free_end] = getticks() - tic_phase;

#endif

  /* Compute group properties and act on the results
   * (seed BHs, dump catalogues..) */
  fof_compute_group_props(e->fof_properties, e->black_holes_properties,
                          e->physical_constants, e->cosmology, e->s,
                          dump_results, dump_debug_results, seed_black_holes);

  /* Reset flag. */
  e->run_fof = 0;

  /* Flag that a FOF has taken place */
  e->step_props |= engine_step_prop_fof;

  /* ... and find the next FOF time */
  if (seed_black_holes) engine_compute_next_fof_time(e);

  /* Free everything we allocated */
  tic_phase = getticks();
  fof_free_arrays(e->fof_properties);
  fof_timers[fof_timer_free_arrays] = getticks() - tic_phase;

  tic_phase = getticks();
#ifdef WITH_MPI
  MPI_Barrier(MPI_COMM_WORLD);
#endif
  fof_timers[fof_timer_end_wait] = getticks() - tic_phase;

  /* Restore the foreign buffers as they were*/
  tic_phase = getticks();
  if (foreign_buffers_allocated) {
#ifdef WITH_MPI
    engine_allocate_foreign_particles(e, /*fof=*/0);
#endif
  }
  fof_timers[fof_timer_realloc_foreign] = getticks() - tic_phase;

  /* TEMPORARY: report the FOF scaling timers (collective). */
  fof_timers[fof_timer_total] = getticks() - tic;
  fof_timers_report(e);

  if (engine_rank == 0)
    message("Complete FOF search took: %.3f %s.",
            clocks_from_ticks(getticks() - tic), clocks_getunit());
#else
  error("SWIFT was not compiled with FOF enabled!");
#endif
}
