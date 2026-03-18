/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2016 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2026 Will J. Roper (w.roper@sussex.ac.uk)
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

/* Some standard headers. */

/* This object's header. */
#include "scheduler.h"

/* Local headers. */
#include "engine.h"
#include "error.h"
#include "sort_part.h"
#include "space.h"
#include "space_getsid.h"
#include "threadpool.h"

/**
 * @brief Split a hydrodynamic task if too large.
 *
 * @param t The #task
 * @param s The #scheduler we are working in.
 */
static void scheduler_splittask_hydro(struct task *t, struct scheduler *s) {
  /* Are we considering both stars and hydro when splitting? */
  /* Note this is not very clean as the scheduler should not really
     access the engine... */
  const int with_feedback = (s->space->e->policy & engine_policy_feedback);
  const int with_stars = (s->space->e->policy & engine_policy_stars);
  const int with_sinks = (s->space->e->policy & engine_policy_sinks);
  const int with_black_holes =
      (s->space->e->policy & engine_policy_black_holes);

  /* Iterate on this task until we're done with it. */
  int redo = 1;
  while (redo) {
    /* Reset the redo flag. */
    redo = 0;

    /* Is this a non-empty self-task? */
    const int is_self =
        (t->type == task_type_self) && (t->ci != NULL) &&
        ((t->ci->hydro.count > 0) || (with_stars && t->ci->stars.count > 0) ||
         (with_sinks && t->ci->sinks.count > 0) ||
         (with_black_holes && t->ci->black_holes.count > 0));

    /* Is this a non-empty pair-task? */
    const int is_pair = (t->type == task_type_pair) && (t->ci != NULL) &&
                        (t->cj != NULL) &&
                        ((t->ci->hydro.count > 0) ||
                         (with_feedback && t->ci->stars.count > 0) ||
                         (with_sinks && t->ci->sinks.count > 0) ||
                         (with_black_holes && t->ci->black_holes.count > 0)) &&
                        ((t->cj->hydro.count > 0) ||
                         (with_feedback && t->cj->stars.count > 0) ||
                         (with_sinks && t->cj->sinks.count > 0) ||
                         (with_black_holes && t->cj->black_holes.count > 0));

    /* Empty task? */
    if (!is_self && !is_pair) {
      t->type = task_type_none;
      t->subtype = task_subtype_none;
      t->ci = NULL;
      t->cj = NULL;
      t->skip = 1;
      break;
    }

    /* Self-interaction? */
    if (t->type == task_type_self) {
      /* Get a handle on the cell involved. */
      struct cell *ci = t->ci;

      /* Foreign task? */
      if (ci->nodeID != s->nodeID) {
        t->skip = 1;
        break;
      }

      /* Is this cell even split and the task does not violate h ? */
      if (cell_can_split_self_hydro_task(ci)) {
        /* Make a sub? */
        if (scheduler_dosub && (ci->hydro.count < space_subsize_self_hydro) &&
            (ci->stars.count < space_subsize_self_stars)) {

          /* Nothing to do here */

          /* Otherwise, make tasks explicitly. */
        } else {
          /* Take a step back (we're going to recycle the current task)... */
          redo = 1;

          /* Add the self tasks. */
          int first_child = 0;
          while (ci->progeny[first_child] == NULL) first_child++;

          t->ci = ci->progeny[first_child];
          cell_set_flag(t->ci, cell_flag_has_tasks);

          for (int k = first_child + 1; k < 8; k++) {
            /* Do we have a non-empty progenitor? */
            if (ci->progeny[k] != NULL &&
                (ci->progeny[k]->hydro.count ||
                 (with_stars && ci->progeny[k]->stars.count))) {
              scheduler_splittask_hydro(
                  scheduler_addtask(s, task_type_self, t->subtype, 0, 0,
                                    ci->progeny[k], NULL),
                  s);
            }
          }

          /* Make a task for each pair of progeny */
          for (int j = 0; j < 8; j++) {
            /* Do we have a non-empty progenitor? */
            if (ci->progeny[j] != NULL &&
                (ci->progeny[j]->hydro.count ||
                 (with_feedback && ci->progeny[j]->stars.count))) {
              for (int k = j + 1; k < 8; k++) {
                /* Do we have a second non-empty progenitor? */
                if (ci->progeny[k] != NULL &&
                    (ci->progeny[k]->hydro.count ||
                     (with_feedback && ci->progeny[k]->stars.count))) {
                  scheduler_splittask_hydro(
                      scheduler_addtask(s, task_type_pair, t->subtype,
                                        sub_sid_flag[j][k], 0, ci->progeny[j],
                                        ci->progeny[k]),
                      s);
                }
              }
            }
          }
        }
      }
    } /* Self interaction */

    /* Pair interaction? */
    else if (t->type == task_type_pair) {
      /* Get a handle on the cells involved. */
      struct cell *ci = t->ci;
      struct cell *cj = t->cj;

      /* Foreign task? */
      if (ci->nodeID != s->nodeID && cj->nodeID != s->nodeID) {
        t->skip = 1;
        break;
      }

      /* Get the sort ID, use space_getsid_and_swap_cells and not t->flags
         to make sure we get ci and cj swapped if needed. */
      double shift[3];
      const int sid = space_getsid_and_swap_cells(s->space, &ci, &cj, shift);

#ifdef SWIFT_DEBUG_CHECKS
      if (sid != t->flags)
        error("Got pair task with incorrect flags: sid=%d flags=%lld", sid,
              t->flags);
#endif

      /* Should this task be split-up? */
      if (cell_can_split_pair_hydro_task(ci) &&
          cell_can_split_pair_hydro_task(cj)) {

        const int h_count_i = ci->hydro.count;
        const int h_count_j = cj->hydro.count;

        const int s_count_i = ci->stars.count;
        const int s_count_j = cj->stars.count;

        int do_sub_hydro = 1;
        int do_sub_stars_i = 1;
        int do_sub_stars_j = 1;
        if (h_count_i > 0 && h_count_j > 0) {

          /* Note: Use division to avoid integer overflow. */
          do_sub_hydro =
              h_count_i * sid_scale[sid] < space_subsize_pair_hydro / h_count_j;
        }
        if (s_count_i > 0 && h_count_j > 0) {

          /* Note: Use division to avoid integer overflow. */
          do_sub_stars_i =
              s_count_i * sid_scale[sid] < space_subsize_pair_stars / h_count_j;
        }
        if (s_count_j > 0 && h_count_i > 0) {

          /* Note: Use division to avoid integer overflow. */
          do_sub_stars_j =
              s_count_j * sid_scale[sid] < space_subsize_pair_stars / h_count_i;
        }

        /* Replace by a single sub-task? */
        if (scheduler_dosub &&
            (do_sub_hydro && do_sub_stars_i && do_sub_stars_j) &&
            !sort_is_corner(sid)) {

          /* Nothing to do here! */

          /* Otherwise, split it. */
        } else {
          /* Take a step back (we're going to recycle the current task)... */
          redo = 1;

          /* Loop over the sub-cell pairs for the current sid and add new tasks
           * for them. */
          struct cell_split_pair *csp = &cell_split_pairs[sid];

          t->ci = ci->progeny[csp->pairs[0].pid];
          t->cj = cj->progeny[csp->pairs[0].pjd];
          if (t->ci != NULL) cell_set_flag(t->ci, cell_flag_has_tasks);
          if (t->cj != NULL) cell_set_flag(t->cj, cell_flag_has_tasks);

          t->flags = csp->pairs[0].sid;
          for (int k = 1; k < csp->count; k++) {
            scheduler_splittask_hydro(
                scheduler_addtask(s, task_type_pair, t->subtype,
                                  csp->pairs[k].sid, 0,
                                  ci->progeny[csp->pairs[k].pid],
                                  cj->progeny[csp->pairs[k].pjd]),
                s);
          }
        }

        /* Otherwise, break it up if it is too large? */
      } else if (scheduler_doforcesplit && ci->split && cj->split &&
                 (ci->hydro.count > space_maxsize / cj->hydro.count)) {
        // message( "force splitting pair with %i and %i parts." ,
        // ci->hydro.count , cj->hydro.count );

        /* Replace the current task. */
        t->type = task_type_none;

        for (int j = 0; j < 8; j++)
          if (ci->progeny[j] != NULL && ci->progeny[j]->hydro.count)
            for (int k = 0; k < 8; k++)
              if (cj->progeny[k] != NULL && cj->progeny[k]->hydro.count) {
                struct task *tl =
                    scheduler_addtask(s, task_type_pair, t->subtype, 0, 0,
                                      ci->progeny[j], cj->progeny[k]);
                scheduler_splittask_hydro(tl, s);
                tl->flags = space_getsid_and_swap_cells(s->space, &t->ci,
                                                        &t->cj, shift);
              }
      }
    } /* pair interaction? */
  } /* iterate over the current task. */
}

/**
 * @brief Split a gravity task if too large.
 *
 * @param t The #task
 * @param s The #scheduler we are working in.
 */
static void scheduler_splittask_gravity(struct task *t, struct scheduler *s) {
  const struct space *sp = s->space;
  struct engine *e = sp->e;

#ifdef SWIFT_DEBUG_CHECKS
  /* Ensure we haven't got a task including a void cell. */
  if (t->ci->subtype == cell_subtype_void ||
      (t->type == task_type_pair && t->cj->subtype == cell_subtype_void))
    error("Got a task with a void cell.");
#endif

#ifdef SWIFT_DEBUG_CHECKS
  /* Mark that we reached this cell in the task splitting. We need to only do
   * this for ci since we check this at the top level (for zoom cells) which
   * should all make it this far for a self task. */
  t->ci->reached_in_task_split = 1;
#endif

  /* Iterate on this task until we're done with it. */
  int redo = 1;
  while (redo) {
    /* Reset the redo flag. */
    redo = 0;

    /* Non-splittable task? */
    if ((t->ci == NULL) || (t->type == task_type_pair && t->cj == NULL)) {
      t->type = task_type_none;
      t->subtype = task_subtype_none;
      t->ci = NULL;
      t->cj = NULL;
      t->skip = 1;
      break;
    }

    /* Self-interaction? */
    if (t->type == task_type_self) {
      /* Get a handle on the cell involved. */
      const struct cell *ci = t->ci;

      /* Foreign task? */
      if (ci->nodeID != s->nodeID) {
        t->skip = 1;
        break;
      }

      /* Should we split this task? */
      if (cell_can_split_self_gravity_task(ci)) {
        if (scheduler_dosub && ci->grav.count < space_subsize_self_grav) {
          /* Otherwise, split it. */
        } else {
          /* Take a step back (we're going to recycle the current task)... */
          redo = 1;

          /* Add the self tasks. */
          int first_child = 0;
          while (ci->progeny[first_child] == NULL) first_child++;

          t->ci = ci->progeny[first_child];
          cell_set_flag(t->ci, cell_flag_has_tasks);

          for (int k = first_child + 1; k < 8; k++)
            if (ci->progeny[k] != NULL)
              scheduler_splittask_gravity(
                  scheduler_addtask(s, task_type_self, t->subtype, 0, 0,
                                    ci->progeny[k], NULL),
                  s);

          /* Make a task for each pair of progeny */
          if (t->subtype != task_subtype_external_grav) {
            for (int j = 0; j < 8; j++) {
              if (ci->progeny[j] == NULL) continue;
              struct cell *cpi = ci->progeny[j];
              for (int k = j + 1; k < 8; k++) {
                if (ci->progeny[k] == NULL) continue;
                struct cell *cpj = ci->progeny[k];

                /* If running with the mesh this pair may be beyond the mesh
                 * criterion meaning we won't need a task here. */
                if (cell_can_use_mesh(e, cpi, cpj)) {
                  continue;
                }

                scheduler_splittask_gravity(
                    scheduler_addtask(s, task_type_pair, t->subtype,
                                      sub_sid_flag[j][k], 0, cpi, cpj),
                    s);
              } /* cpj loop */
            } /* cpi loop */
          } /* Self-gravity only */
        } /* Make tasks explicitly */
      } /* Cell is split */
    } /* Self interaction */

    /* Pair interaction? */
    else if (t->type == task_type_pair) {
      /* Get a handle on the cells involved. */
      struct cell *ci = t->ci;
      struct cell *cj = t->cj;

      /* Foreign task? */
      if (ci->nodeID != s->nodeID && cj->nodeID != s->nodeID) {
        t->skip = 1;
        break;
      }

      /* Should this task be split-up? */
      if (cell_can_split_pair_gravity_task(ci, cj)) {
        const long long gcount_i = ci->grav.count;
        const long long gcount_j = cj->grav.count;

        /* Replace by a single sub-task? */
        if (scheduler_dosub &&
            gcount_i * gcount_j < ((long long)space_subsize_pair_grav)) {
          /* Otherwise, split it. */
        } else {
          /* Turn the task into a M-M task that will take care of all the
           * progeny pairs */
          t->type = task_type_grav_mm;
          t->subtype = task_subtype_progeny;
          t->flags = 0;

          /* Make a task for every other pair of progeny */
          for (int i = 0; i < 8; i++) {
            if (ci->progeny[i] != NULL) {
              struct cell *cpi = ci->progeny[i];
              for (int j = 0; j < 8; j++) {
                if (cj->progeny[j] != NULL) {
                  struct cell *cpj = cj->progeny[j];

                  /* If running with the mesh this pair may be beyond the mesh
                   * criterion meaning we won't need a task here. */
                  if (cell_can_use_mesh(e, cpi, cpj)) {
                    continue;
                  }

                  /* Can we use a M-M interaction here? */
                  if (cell_can_use_pair_mm(ci->progeny[i], cj->progeny[j], e,
                                           sp, /*use_rebuild_data=*/1,
                                           /*is_tree_walk=*/1,
                                           /*periodic boundaries*/ sp->periodic,
                                           /*use_mesh*/ sp->periodic)) {

                    /* Flag this pair as being treated by the M-M task.
                     * We use the 64 bits in the task->flags field to store
                     * this information. The corresponding task will unpack
                     * the information and operate according to the choices
                     * made here. */
                    const int flag = i * 8 + j;
                    t->flags |= (1ULL << flag);

                  } else {
                    /* Ok, we actually have to create a task */
                    scheduler_splittask_gravity(
                        scheduler_addtask(s, task_type_pair, task_subtype_grav,
                                          0, 0, cpi, cpj),
                        s);
                  }
                }
              }
            }
          }

          /* Can none of the progenies use M-M calculations? */
          if (t->flags == 0) {
            t->type = task_type_none;
            t->subtype = task_subtype_none;
            t->ci = NULL;
            t->cj = NULL;
            t->skip = 1;
          }

        } /* Split the pair */
      }
    } /* pair interaction? */
  } /* iterate over the current task. */
}

/**
 * @brief Split a void cell pair gravity task to get down to the zoom cells.
 *
 * Both cells must be either split or void (void cells are always split but
 * have split=0 at the zoom interface). If neither cell is void the task is
 * redirected to the normal gravity splitter. Otherwise the task is converted
 * to a grav_mm/progeny task that handles M-M interactions between all pairs
 * of progeny, with any progeny pairs that cannot use M-M spawned as new pair
 * tasks and recursed into. Once the zoom region is reached
 * scheduler_splittask_gravity handles any further splitting.
 *
 * @param t The #task
 * @param s The #scheduler we are working in.
 */
static void zoom_scheduler_splittask_gravity_void_pair(struct task *t,
                                                       struct scheduler *s) {
  const struct space *sp = s->space;
  struct engine *e = sp->e;

#ifdef SWIFT_DEBUG_CHECKS
  /* Ensure we have a pair task. */
  if (t->type != task_type_pair) {
    error("Got a non-pair task (t->type=%s)", taskID_names[t->type]);
  }

  /* Ensure everything is at the same depth (while this is the norm, we
   * once didn't split symmetrically so should make sure this is
   * not the case here). */
  int depth_i = t->ci->depth;
  int depth_j = t->cj->depth;

  if (t->ci->type == cell_type_zoom) {
    depth_i += sp->zoom_props->zoom_cell_depth;
  }
  if (t->cj->type == cell_type_zoom) {
    depth_j += sp->zoom_props->zoom_cell_depth;
  }
  if (depth_i != depth_j) {
    error("Got a pair task with different depths: %d (%s/%s) != %d (%s/%s)",
          depth_i, cellID_names[t->ci->type], subcellID_names[t->ci->subtype],
          depth_j, cellID_names[t->cj->type], subcellID_names[t->cj->subtype]);
  }

#endif

  /* Get a handle on the cells involved. */
  struct cell *restrict ci = t->ci;
  struct cell *restrict cj = t->cj;

  /* If neither cell is a void cell, redirect to the normal splitter. */
  if (ci->subtype != cell_subtype_void && cj->subtype != cell_subtype_void) {
    scheduler_splittask_gravity(t, s);
    return;
  }

  /* Both cells must be splittable: either genuinely split, or a void cell
   * (which always has progeny despite having split=0 at the zoom interface).
   * If either side is an unsplit non-void leaf we cannot iterate progeny. */
  if (!(ci->split || ci->subtype == cell_subtype_void) ||
      !(cj->split || cj->subtype == cell_subtype_void)) {
    error(
        "zoom_scheduler_splittask_gravity_void_pair called with a non-split "
        "non-void cell: ci->split=%d ci->subtype=%s cj->split=%d "
        "cj->subtype=%s",
        ci->split, subcellID_names[ci->subtype], cj->split,
        subcellID_names[cj->subtype]);
  }

  /* Convert to a grav_mm/progeny task. The flags field encodes which progeny
   * pairs will be handled by M-M; remaining pairs are spawned as new tasks. */
  t->type = task_type_grav_mm;
  t->subtype = task_subtype_progeny;
  t->flags = 0;

  /* Define a flag for when the original task has been reused. */
  int reused = 0;

  /* When we split a regular cell's task because it is interacting with a
   * void cell, we can end up below the depth set by space_subdepth_diff_grav.
   * This will cause absolute havoc with hierarchical gravity tasks being
   * missing on the regular cell if we don't flag this somehow to ensure
   * task recursions continue to this level. */
  if (!cell_is_above_diff_grav_depth(ci)) {
    ci->grav.tasks_below_diff_grav_depth = 1;
  }
  if (!cell_is_above_diff_grav_depth(cj)) {
    cj->grav.tasks_below_diff_grav_depth = 1;
  }

  /* Loop over the progeny. */
  for (int i = 0; i < 8; i++) {
    struct cell *cpi = ci->progeny[i];

    /* Skip NULL progeny (can happen if a background cell tree has no
     * particles at the leaves for whatever reason). */
    if (cpi == NULL) continue;

    /* Skip any empty progeny of a void cell (void cells themselves always
     * have 0 particles but are never "empty"). */
    if (cell_is_empty_grav(cpi)) continue;

    for (int j = 0; j < 8; j++) {
      struct cell *cpj = cj->progeny[j];

      /* Skip NULL progeny (can happen is a background cell tree has no
       * particles at the leaves for whatever reason). */
      if (cpj == NULL) continue;

      /* Skip any empty progeny of a void cell (void cells themselves always
       * have 0 particles but are never "empty"). */
      if (cell_is_empty_grav(cpj)) continue;

      /* Skip entirely foreign pairs. */
      if (cpi->nodeID != engine_rank && cpj->nodeID != engine_rank) continue;

      /* Could we use the mesh for this pair? */
      if (cell_can_use_mesh(e, cpi, cpj)) {
        continue;
      }

      /* Can we use a M-M interaction here? */
      if (cell_can_use_pair_mm(cpi, cpj, e, sp,
                               /*use_rebuild_data=*/1,
                               /*is_tree_walk=*/1,
                               /*periodic boundaries*/ sp->periodic,
                               /*use_mesh*/ sp->periodic)) {

        /* Flag that we've reused the original task. */
        reused = 1;

        /* Flag this pair as being treated by the M-M task.
         * We use the 64 bits in the task->flags field to store
         * this information. The corresponding taks will unpack
         * the information and operate according to the choices
         * made here. */
        const int flag = i * 8 + j;
        t->flags |= (1ULL << flag);

      } else {

        /* Can't use an M-M so let's make a pair task. */
        zoom_scheduler_splittask_gravity_void_pair(
            scheduler_addtask(s, task_type_pair, task_subtype_grav, 0, 0, cpi,
                              cpj),
            s);
      }
    }
  }

  /* Can none of the progenies use M-M calculations? */
  if (!reused) {
    t->type = task_type_none;
    t->subtype = task_subtype_none;
    t->ci = NULL;
    t->cj = NULL;
    t->skip = 1;
  }
}

/**
 * @brief Split a void cell self gravity task to get down to the zoom cells.
 *
 * This will create a self task for each progeny and a pair task for each
 * progeny pair. This will continue splitting until we reach the zoom cells,
 * at which point scheduler_splittask_gravity will take over.
 *
 * @param t The #task
 * @param s The #scheduler we are working in.
 */
static void zoom_scheduler_splittask_gravity_void_self(struct task *t,
                                                       struct scheduler *s) {

#ifdef SWIFT_DEBUG_CHECKS
  /* Ensure we have a self task. */
  if (t->type != task_type_self) {
    error("Got a non-self task (t->type=%s)", taskID_names[t->type]);
  }
#endif

  /* Iterate on this task until we're done with it. */
  while (t->ci->subtype == cell_subtype_void) {

    /* Get a handle on the cell involved. */
    const struct cell *ci = t->ci;

    /* Get the first progeny that exists and is local (and in the case
     * of a non-void progeny: has particles. */
    int first_child = 0;
    for (; first_child < 8; first_child++) {
      /* Skip non-existant progeny. */
      if (ci->progeny[first_child] == NULL) continue;

      /* Skip empty progeny. */
      if (cell_is_empty_grav(ci->progeny[first_child])) continue;

      /* Skip foreign zoom progeny (no such thing as a foreign self task). */
      if (ci->progeny[first_child]->type == cell_type_zoom &&
          ci->progeny[first_child]->nodeID != engine_rank) {
        continue;
      }

      /* We found our first child. Stop looking. */
      break;
    }

    /* If we have found no progeny, we are done here. */
    if (first_child == 8) break;

    /* Reuse the task we already have. */
    t->ci = ci->progeny[first_child];
    cell_set_flag(t->ci, cell_flag_has_tasks);

    /* Create a self for all progeny beyond the first. */
    for (int i = first_child + 1; i < 8; i++) {

      /* Skip non-existant progeny. */
      if (ci->progeny[i] == NULL) continue;

      /* Skip empty progeny. */
      if (cell_is_empty_grav(ci->progeny[i])) continue;

      /* Skip non-local progeny (no such thing as a foreign self task). */
      if (ci->progeny[i]->type == cell_type_zoom &&
          ci->progeny[i]->nodeID != engine_rank)
        continue;

      /* Create the self task. */
      zoom_scheduler_splittask_gravity_void_self(
          scheduler_addtask(s, task_type_self, t->subtype, 0, 0, ci->progeny[i],
                            NULL),
          s);
    }

    /* Create pair tasks for all pairs of progeny. */
    for (int j = 0; j < 8; j++) {

      /* Skip non-existant progeny. */
      if (ci->progeny[j] == NULL) continue;

      /* Skip empty non-void progeny. */
      if (ci->progeny[j]->subtype != cell_subtype_void &&
          ci->progeny[j]->grav.count == 0)
        continue;

      for (int k = j + 1; k < 8; k++) {

        /* Skip non-existant progeny. */
        if (ci->progeny[k] == NULL) continue;

        /* Skip empty progeny. */
        if (cell_is_empty_grav(ci->progeny[k])) continue;

        /* Skip entirely foreign pairs. */
        if ((ci->progeny[j]->type == cell_type_zoom &&
             ci->progeny[k]->type == cell_type_zoom) &&
            ci->progeny[j]->nodeID != engine_rank &&
            ci->progeny[k]->nodeID != engine_rank)
          continue;

        /* Could we use the mesh for this pair? */
        if (cell_can_use_mesh(s->space->e, ci->progeny[j], ci->progeny[k])) {
          continue;
        }

        /* Create the pair task. */
        zoom_scheduler_splittask_gravity_void_pair(
            scheduler_addtask(s, task_type_pair, t->subtype, sub_sid_flag[j][k],
                              0, ci->progeny[j], ci->progeny[k]),
            s);
      }
    }
  }

  /* Exit and kill the task if we ended up with an empty task, a foreign task or
   * we still have a void cell. */
  if (t->ci == NULL || t->ci->nodeID != engine_rank ||
      t->ci->subtype == cell_subtype_void) {
    t->type = task_type_none;
    t->subtype = task_subtype_none;
    t->ci = NULL;
    t->skip = 1;
    return;
  }

  /* Now we're not in a void cell we can just call the normal splitter.  */
  scheduler_splittask_gravity(t, s);
}

/**
 * @brief Split a FOF task if too large.
 *
 * @param t The #task
 * @param s The #scheduler we are working in.
 */
static void scheduler_splittask_fof(struct task *t, struct scheduler *s) {

  /* Iterate on this task until we're done with it. */
  int redo = 1;
  while (redo) {

    /* Reset the redo flag. */
    redo = 0;

    /* Non-splittable task? */
    if ((t->ci == NULL) || (t->type == task_type_fof_pair && t->cj == NULL) ||
        t->ci->grav.count == 0 || (t->cj != NULL && t->cj->grav.count == 0)) {
      t->type = task_type_none;
      t->subtype = task_subtype_none;
      t->ci = NULL;
      t->cj = NULL;
      t->skip = 1;
      break;
    }

    /* Self-interaction? */
    if (t->type == task_type_fof_self) {

      /* Get a handle on the cell involved. */
      struct cell *ci = t->ci;

      /* Foreign task? */
      if (ci->nodeID != s->nodeID) {
        t->skip = 1;
        break;
      }

      /* Is this cell even split? */
      if (cell_can_split_self_fof_task(ci)) {

        /* Take a step back (we're going to recycle the current task)... */
        redo = 1;

        /* Add the self tasks. */
        int first_child = 0;
        while (ci->progeny[first_child] == NULL) first_child++;
        t->ci = ci->progeny[first_child];
        for (int k = first_child + 1; k < 8; k++)
          if (ci->progeny[k] != NULL && ci->progeny[k]->grav.count)
            scheduler_splittask_fof(
                scheduler_addtask(s, task_type_fof_self, t->subtype, 0, 0,
                                  ci->progeny[k], NULL),
                s);

        /* Make a task for each pair of progeny */
        for (int j = 0; j < 8; j++)
          if (ci->progeny[j] != NULL && ci->progeny[j]->grav.count)
            for (int k = j + 1; k < 8; k++)
              if (ci->progeny[k] != NULL && ci->progeny[k]->grav.count)
                scheduler_splittask_fof(
                    scheduler_addtask(s, task_type_fof_pair, t->subtype, 0, 0,
                                      ci->progeny[j], ci->progeny[k]),
                    s);
      } /* Cell is split */

    } /* Self interaction */

  } /* iterate over the current task. */
}

/**
 * @brief Mapper function to split FOF tasks that may be too large.
 *
 * @param map_data the tasks to process
 * @param num_elements the number of tasks.
 * @param extra_data The #scheduler we are working in.
 */
void scheduler_splittasks_fof_mapper(void *map_data, int num_elements,
                                     void *extra_data) {
  /* Extract the parameters. */
  struct scheduler *s = (struct scheduler *)extra_data;
  struct task *tasks = (struct task *)map_data;

  for (int ind = 0; ind < num_elements; ind++) {
    struct task *t = &tasks[ind];

    /* Invoke the correct splitting strategy */
    if (t->type == task_type_fof_self || t->type == task_type_fof_pair) {
      scheduler_splittask_fof(t, s);
    }
  }
}

/**
 * @brief Mapper function to split non-FOF tasks that may be too large.
 *
 * @param map_data the tasks to process
 * @param num_elements the number of tasks.
 * @param extra_data The #scheduler we are working in.
 */
void scheduler_splittasks_mapper(void *map_data, int num_elements,
                                 void *extra_data) {
  /* Extract the parameters. */
  struct scheduler *s = (struct scheduler *)extra_data;
  struct task *tasks = (struct task *)map_data;

  for (int ind = 0; ind < num_elements; ind++) {
    struct task *t = &tasks[ind];

    /* Invoke the correct splitting strategy */
    if (t->subtype == task_subtype_density) {
      scheduler_splittask_hydro(t, s);
    } else if (t->subtype == task_subtype_external_grav) {
      scheduler_splittask_gravity(t, s);
    } else if (t->type == task_type_self &&
               t->ci->subtype == cell_subtype_void) {
      zoom_scheduler_splittask_gravity_void_self(t, s);
    } else if (t->type == task_type_pair &&
               (t->ci->subtype == cell_subtype_void ||
                t->cj->subtype == cell_subtype_void)) {
      zoom_scheduler_splittask_gravity_void_pair(t, s);
    } else if (t->subtype == task_subtype_grav) {
      scheduler_splittask_gravity(t, s);
    } else {
#ifdef SWIFT_DEBUG_CHECKS
      error("Unexpected task sub-type %s/%s", taskID_names[t->type],
            subtaskID_names[t->subtype]);
#endif
    }
  }
}

/**
 * @brief Splits all the tasks in the scheduler that are too large.
 *
 * @param s The #scheduler.
 * @param fof_tasks Are we splitting the FOF tasks (1)? Or the regular tasks
 * (0)?
 * @param verbose Are we talkative?
 */
void scheduler_splittasks(struct scheduler *s, const int fof_tasks,
                          const int verbose) {

  if (verbose) {
    message("space_subsize_self_hydro= %d", space_subsize_self_hydro);
    message("space_subsize_pair_hydro= %d", space_subsize_pair_hydro);
    message("space_subsize_self_stars= %d", space_subsize_self_stars);
    message("space_subsize_pair_stars= %d", space_subsize_pair_stars);
    message("space_subsize_self_grav= %d", space_subsize_self_grav);
    message("space_subsize_pair_grav= %d", space_subsize_pair_grav);
  }

  if (fof_tasks) {
    /* Call the mapper on each current task. */
    threadpool_map(s->threadpool, scheduler_splittasks_fof_mapper, s->tasks,
                   s->nr_tasks, sizeof(struct task), threadpool_auto_chunk_size,
                   s);

  } else {
    /* Call the mapper on each current task. */
    threadpool_map(s->threadpool, scheduler_splittasks_mapper, s->tasks,
                   s->nr_tasks, sizeof(struct task), threadpool_auto_chunk_size,
                   s);
  }
}
