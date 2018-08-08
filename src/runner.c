/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
 *               2016 John A. Regan (john.a.regan@durham.ac.uk)
 *                    Tom Theuns (tom.theuns@durham.ac.uk)
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

/* Some standard headers. */
#include <float.h>
#include <limits.h>
#include <stdlib.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "runner.h"

/* Local headers. */
#include "active.h"
#include "approx_math.h"
#include "atomic.h"
#include "cell.h"
#include "chemistry.h"
#include "const.h"
#include "cooling.h"
#include "debug.h"
#include "drift.h"
#include "engine.h"
#include "error.h"
#include "gravity.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "kick.h"
#include "minmax.h"
#include "runner_doiact_vec.h"
#include "scheduler.h"
#include "sort_part.h"
#include "sourceterms.h"
#include "space.h"
#include "space_getsid.h"
#include "stars.h"
#include "task.h"
#include "timers.h"
#include "timestep.h"

#define TASK_LOOP_DENSITY 0
#define TASK_LOOP_GRADIENT 1
#define TASK_LOOP_FORCE 2
#define TASK_LOOP_LIMITER 3

/* Import the density loop functions. */
#define FUNCTION density
#define FUNCTION_TASK_LOOP TASK_LOOP_DENSITY
#include "runner_doiact.h"
#undef FUNCTION
#undef FUNCTION_TASK_LOOP

/* Import the gradient loop functions (if required). */
#ifdef EXTRA_HYDRO_LOOP
#define FUNCTION gradient
#define FUNCTION_TASK_LOOP TASK_LOOP_GRADIENT
#include "runner_doiact.h"
#undef FUNCTION
#undef FUNCTION_TASK_LOOP
#endif

/* Import the force loop functions. */
#define FUNCTION force
#define FUNCTION_TASK_LOOP TASK_LOOP_FORCE
#include "runner_doiact.h"
#undef FUNCTION
#undef FUNCTION_TASK_LOOP

/* Import the gravity loop functions. */
#include "runner_doiact_grav.h"

/**
 * @brief Perform source terms
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void runner_do_sourceterms(struct runner *r, struct cell *c, int timer) {
  const int count = c->count;
  const double cell_min[3] = {c->loc[0], c->loc[1], c->loc[2]};
  const double cell_width[3] = {c->width[0], c->width[1], c->width[2]};
  struct sourceterms *sourceterms = r->e->sourceterms;
  const int dimen = 3;

  TIMER_TIC;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_sourceterms(r, c->progeny[k], 0);
  } else {

    if (count > 0) {

      /* do sourceterms in this cell? */
      const int incell =
          sourceterms_test_cell(cell_min, cell_width, sourceterms, dimen);
      if (incell == 1) {
        sourceterms_apply(r, sourceterms, c);
      }
    }
  }

  if (timer) TIMER_TOC(timer_dosource);
}

/**
 * @brief Calculate gravity acceleration from external potential
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void runner_do_grav_external(struct runner *r, struct cell *c, int timer) {

  struct gpart *restrict gparts = c->gparts;
  const int gcount = c->gcount;
  const struct engine *e = r->e;
  const struct external_potential *potential = e->external_potential;
  const struct phys_const *constants = e->physical_constants;
  const double time = r->e->time;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_gravity(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_grav_external(r, c->progeny[k], 0);
  } else {

    /* Loop over the gparts in this cell. */
    for (int i = 0; i < gcount; i++) {

      /* Get a direct pointer on the part. */
      struct gpart *restrict gp = &gparts[i];

      /* Is this part within the time step? */
      if (gpart_is_active(gp, e)) {
        external_gravity_acceleration(time, potential, constants, gp);
      }
    }
  }

  if (timer) TIMER_TOC(timer_dograv_external);
}

/**
 * @brief Calculate gravity accelerations from the periodic mesh
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void runner_do_grav_mesh(struct runner *r, struct cell *c, int timer) {

  struct gpart *restrict gparts = c->gparts;
  const int gcount = c->gcount;
  const struct engine *e = r->e;

#ifdef SWIFT_DEBUG_CHECKS
  if (!e->s->periodic) error("Calling mesh forces in non-periodic mode.");
#endif

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_gravity(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_grav_mesh(r, c->progeny[k], 0);
  } else {

    /* Get the forces from the gravity mesh */
    pm_mesh_interpolate_forces(e->mesh, e, gparts, gcount);
  }

  if (timer) TIMER_TOC(timer_dograv_mesh);
}

/**
 * @brief Calculate change in thermal state of particles induced
 * by radiative cooling and heating.
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void runner_do_cooling(struct runner *r, struct cell *c, int timer) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const struct cooling_function_data *cooling_func = e->cooling_func;
  const struct phys_const *constants = e->physical_constants;
  const struct unit_system *us = e->internal_units;
  const double time_base = e->time_base;
  const integertime_t ti_current = e->ti_current;
  struct part *restrict parts = c->parts;
  struct xpart *restrict xparts = c->xparts;
  const int count = c->count;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_hydro(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_cooling(r, c->progeny[k], 0);
  } else {

    /* Loop over the parts in this cell. */
    for (int i = 0; i < count; i++) {

      /* Get a direct pointer on the part. */
      struct part *restrict p = &parts[i];
      struct xpart *restrict xp = &xparts[i];

      if (part_is_active(p, e)) {

        double dt_cool;
        if (with_cosmology) {
          const integertime_t ti_step = get_integer_timestep(p->time_bin);
          const integertime_t ti_begin =
              get_integer_time_begin(ti_current + 1, p->time_bin);
          dt_cool =
              cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
        } else {
          dt_cool = get_timestep(p->time_bin, time_base);
        }

        /* Let's cool ! */
        cooling_cool_part(constants, us, cosmo, cooling_func, p, xp, dt_cool);
      }
    }
  }

  if (timer) TIMER_TOC(timer_do_cooling);
}

/**
 * @brief Sort the entries in ascending order using QuickSort.
 *
 * @param sort The entries
 * @param N The number of entries.
 */
void runner_do_sort_ascending(struct entry *sort, int N) {

  struct {
    short int lo, hi;
  } qstack[10];
  int qpos, i, j, lo, hi, imin;
  struct entry temp;
  float pivot;

  /* Sort parts in cell_i in decreasing order with quicksort */
  qstack[0].lo = 0;
  qstack[0].hi = N - 1;
  qpos = 0;
  while (qpos >= 0) {
    lo = qstack[qpos].lo;
    hi = qstack[qpos].hi;
    qpos -= 1;
    if (hi - lo < 15) {
      for (i = lo; i < hi; i++) {
        imin = i;
        for (j = i + 1; j <= hi; j++)
          if (sort[j].d < sort[imin].d) imin = j;
        if (imin != i) {
          temp = sort[imin];
          sort[imin] = sort[i];
          sort[i] = temp;
        }
      }
    } else {
      pivot = sort[(lo + hi) / 2].d;
      i = lo;
      j = hi;
      while (i <= j) {
        while (sort[i].d < pivot) i++;
        while (sort[j].d > pivot) j--;
        if (i <= j) {
          if (i < j) {
            temp = sort[i];
            sort[i] = sort[j];
            sort[j] = temp;
          }
          i += 1;
          j -= 1;
        }
      }
      if (j > (lo + hi) / 2) {
        if (lo < j) {
          qpos += 1;
          qstack[qpos].lo = lo;
          qstack[qpos].hi = j;
        }
        if (i < hi) {
          qpos += 1;
          qstack[qpos].lo = i;
          qstack[qpos].hi = hi;
        }
      } else {
        if (i < hi) {
          qpos += 1;
          qstack[qpos].lo = i;
          qstack[qpos].hi = hi;
        }
        if (lo < j) {
          qpos += 1;
          qstack[qpos].lo = lo;
          qstack[qpos].hi = j;
        }
      }
    }
  }
}

/**
 * @brief Recursively checks that the flags are consistent in a cell hierarchy.
 *
 * Debugging function.
 *
 * @param c The #cell to check.
 * @param flags The sorting flags to check.
 */
void runner_check_sorts(struct cell *c, int flags) {

#ifdef SWIFT_DEBUG_CHECKS
  if (flags & ~c->sorted) error("Inconsistent sort flags (downward)!");
  if (c->split)
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL && c->progeny[k]->count > 0)
        runner_check_sorts(c->progeny[k], c->sorted);
#else
  error("Calling debugging code without debugging flag activated.");
#endif
}

/**
 * @brief Sort the particles in the given cell along all cardinal directions.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param flags Cell flag.
 * @param cleanup If true, re-build the sorts for the selected flags instead
 *        of just adding them.
 * @param clock Flag indicating whether to record the timing or not, needed
 *      for recursive calls.
 */
void runner_do_sort(struct runner *r, struct cell *c, int flags, int cleanup,
                    int clock) {

  struct entry *fingers[8];
  const int count = c->count;
  const struct part *parts = c->parts;
  struct xpart *xparts = c->xparts;
  float buff[8];

  TIMER_TIC;

  /* We need to do the local sorts plus whatever was requested further up. */
  flags |= c->do_sort;
  if (cleanup) {
    c->sorted = 0;
  } else {
    flags &= ~c->sorted;
  }
  if (flags == 0 && !c->do_sub_sort) return;

  /* Check that the particles have been moved to the current time */
  if (flags && !cell_are_part_drifted(c, r->e))
    error("Sorting un-drifted cell c->nodeID=%d", c->nodeID);

#ifdef SWIFT_DEBUG_CHECKS
  /* Make sure the sort flags are consistent (downward). */
  runner_check_sorts(c, c->sorted);

  /* Make sure the sort flags are consistent (upard). */
  for (struct cell *finger = c->parent; finger != NULL;
       finger = finger->parent) {
    if (finger->sorted & ~c->sorted) error("Inconsistent sort flags (upward).");
  }

  /* Update the sort timer which represents the last time the sorts
     were re-set. */
  if (c->sorted == 0) c->ti_sort = r->e->ti_current;
#endif

  /* start by allocating the entry arrays in the requested dimensions. */
  for (int j = 0; j < 13; j++) {
    if ((flags & (1 << j)) && c->sort[j] == NULL) {
      if ((c->sort[j] = (struct entry *)malloc(sizeof(struct entry) *
                                               (count + 1))) == NULL)
        error("Failed to allocate sort memory.");
    }
  }

  /* Does this cell have any progeny? */
  if (c->split) {

    /* Fill in the gaps within the progeny. */
    float dx_max_sort = 0.0f;
    float dx_max_sort_old = 0.0f;
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL && c->progeny[k]->count > 0) {
        /* Only propagate cleanup if the progeny is stale. */
        runner_do_sort(r, c->progeny[k], flags,
                       cleanup && (c->progeny[k]->dx_max_sort >
                                   space_maxreldx * c->progeny[k]->dmin),
                       0);
        dx_max_sort = max(dx_max_sort, c->progeny[k]->dx_max_sort);
        dx_max_sort_old = max(dx_max_sort_old, c->progeny[k]->dx_max_sort_old);
      }
    }
    c->dx_max_sort = dx_max_sort;
    c->dx_max_sort_old = dx_max_sort_old;

    /* Loop over the 13 different sort arrays. */
    for (int j = 0; j < 13; j++) {

      /* Has this sort array been flagged? */
      if (!(flags & (1 << j))) continue;

      /* Init the particle index offsets. */
      int off[8];
      off[0] = 0;
      for (int k = 1; k < 8; k++)
        if (c->progeny[k - 1] != NULL)
          off[k] = off[k - 1] + c->progeny[k - 1]->count;
        else
          off[k] = off[k - 1];

      /* Init the entries and indices. */
      int inds[8];
      for (int k = 0; k < 8; k++) {
        inds[k] = k;
        if (c->progeny[k] != NULL && c->progeny[k]->count > 0) {
          fingers[k] = c->progeny[k]->sort[j];
          buff[k] = fingers[k]->d;
          off[k] = off[k];
        } else
          buff[k] = FLT_MAX;
      }

      /* Sort the buffer. */
      for (int i = 0; i < 7; i++)
        for (int k = i + 1; k < 8; k++)
          if (buff[inds[k]] < buff[inds[i]]) {
            int temp_i = inds[i];
            inds[i] = inds[k];
            inds[k] = temp_i;
          }

      /* For each entry in the new sort list. */
      struct entry *finger = c->sort[j];
      for (int ind = 0; ind < count; ind++) {

        /* Copy the minimum into the new sort array. */
        finger[ind].d = buff[inds[0]];
        finger[ind].i = fingers[inds[0]]->i + off[inds[0]];

        /* Update the buffer. */
        fingers[inds[0]] += 1;
        buff[inds[0]] = fingers[inds[0]]->d;

        /* Find the smallest entry. */
        for (int k = 1; k < 8 && buff[inds[k]] < buff[inds[k - 1]]; k++) {
          int temp_i = inds[k - 1];
          inds[k - 1] = inds[k];
          inds[k] = temp_i;
        }

      } /* Merge. */

      /* Add a sentinel. */
      c->sort[j][count].d = FLT_MAX;
      c->sort[j][count].i = 0;

      /* Mark as sorted. */
      atomic_or(&c->sorted, 1 << j);

    } /* loop over sort arrays. */

  } /* progeny? */

  /* Otherwise, just sort. */
  else {

    /* Reset the sort distance */
    if (c->sorted == 0) {
#ifdef SWIFT_DEBUG_CHECKS
      if (xparts != NULL && c->nodeID != engine_rank)
        error("Have non-NULL xparts in foreign cell");
#endif

      /* And the individual sort distances if we are a local cell */
      if (xparts != NULL) {
        for (int k = 0; k < count; k++) {
          xparts[k].x_diff_sort[0] = 0.0f;
          xparts[k].x_diff_sort[1] = 0.0f;
          xparts[k].x_diff_sort[2] = 0.0f;
        }
      }
      c->dx_max_sort_old = 0.f;
      c->dx_max_sort = 0.f;
    }

    /* Fill the sort array. */
    for (int k = 0; k < count; k++) {
      const double px[3] = {parts[k].x[0], parts[k].x[1], parts[k].x[2]};
      for (int j = 0; j < 13; j++)
        if (flags & (1 << j)) {
          c->sort[j][k].i = k;
          c->sort[j][k].d = px[0] * runner_shift[j][0] +
                            px[1] * runner_shift[j][1] +
                            px[2] * runner_shift[j][2];
        }
    }

    /* Add the sentinel and sort. */
    for (int j = 0; j < 13; j++)
      if (flags & (1 << j)) {
        c->sort[j][count].d = FLT_MAX;
        c->sort[j][count].i = 0;
        runner_do_sort_ascending(c->sort[j], count);
        atomic_or(&c->sorted, 1 << j);
      }
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify the sorting. */
  for (int j = 0; j < 13; j++) {
    if (!(flags & (1 << j))) continue;
    struct entry *finger = c->sort[j];
    for (int k = 1; k < count; k++) {
      if (finger[k].d < finger[k - 1].d)
        error("Sorting failed, ascending array.");
      if (finger[k].i >= count) error("Sorting failed, indices borked.");
    }
  }

  /* Make sure the sort flags are consistent (downward). */
  runner_check_sorts(c, flags);

  /* Make sure the sort flags are consistent (upward). */
  for (struct cell *finger = c->parent; finger != NULL;
       finger = finger->parent) {
    if (finger->sorted & ~c->sorted) error("Inconsistent sort flags.");
  }
#endif

  /* Clear the cell's sort flags. */
  c->do_sort = 0;
  c->do_sub_sort = 0;
  c->requires_sorts = 0;

  if (clock) TIMER_TOC(timer_dosort);
}

/**
 * @brief Initialize the multipoles before the gravity calculation.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer 1 if the time is to be recorded.
 */
void runner_do_init_grav(struct runner *r, struct cell *c, int timer) {

  const struct engine *e = r->e;

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (!(e->policy & engine_policy_self_gravity))
    error("Grav-init task called outside of self-gravity calculation");
#endif

  /* Anything to do here? */
  if (!cell_is_active_gravity(c, e)) return;

  /* Reset the gravity acceleration tensors */
  gravity_field_tensors_init(&c->multipole->pot, e->ti_current);

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) runner_do_init_grav(r, c->progeny[k], 0);
    }
  }

  if (timer) TIMER_TOC(timer_init_grav);
}

/**
 * @brief Intermediate task after the gradient loop that does final operations
 * on the gradient quantities and optionally slope limits the gradients
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_extra_ghost(struct runner *r, struct cell *c, int timer) {

#ifdef EXTRA_HYDRO_LOOP

  struct part *restrict parts = c->parts;
  struct xpart *restrict xparts = c->xparts;
  const int count = c->count;
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_hydro(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_extra_ghost(r, c->progeny[k], 0);
  } else {

    /* Loop over the parts in this cell. */
    for (int i = 0; i < count; i++) {

      /* Get a direct pointer on the part. */
      struct part *restrict p = &parts[i];
      struct xpart *restrict xp = &xparts[i];

      if (part_is_active(p, e)) {

        /* Finish the gradient calculation */
        hydro_end_gradient(p);

        /* As of here, particle force variables will be set. */

        /* Compute variables required for the force loop */
        hydro_prepare_force(p, xp, cosmo);

        /* The particle force values are now set.  Do _NOT_
           try to read any particle density variables! */

        /* Prepare the particle for the force loop over neighbours */
        hydro_reset_acceleration(p);
      }
    }
  }

  if (timer) TIMER_TOC(timer_do_extra_ghost);

#else
  error("SWIFT was not compiled with the extra hydro loop activated.");
#endif
}

/**
 * @brief Intermediate task after the density to check that the smoothing
 * lengths are correct.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_ghost(struct runner *r, struct cell *c, int timer) {

  struct part *restrict parts = c->parts;
  struct xpart *restrict xparts = c->xparts;
  const struct engine *e = r->e;
  const struct space *s = e->s;
  const struct hydro_space *hs = &s->hs;
  const struct cosmology *cosmo = e->cosmology;
  const struct chemistry_global_data *chemistry = e->chemistry;
  const float hydro_h_max = e->hydro_properties->h_max;
  const float eps = e->hydro_properties->h_tolerance;
  const float hydro_eta_dim =
      pow_dimension(e->hydro_properties->eta_neighbours);
  const int max_smoothing_iter = e->hydro_properties->max_smoothing_iterations;
  int redo = 0, count = 0;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_hydro(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_ghost(r, c->progeny[k], 0);
  } else {

    /* Init the list of active particles that have to be updated. */
    int *pid = NULL;
    if ((pid = (int *)malloc(sizeof(int) * c->count)) == NULL)
      error("Can't allocate memory for pid.");
    for (int k = 0; k < c->count; k++)
      if (part_is_active(&parts[k], e)) {
        pid[count] = k;
        ++count;
      }

    /* While there are particles that need to be updated... */
    for (int num_reruns = 0; count > 0 && num_reruns < max_smoothing_iter;
         num_reruns++) {

      /* Reset the redo-count. */
      redo = 0;

      /* Loop over the remaining active parts in this cell. */
      for (int i = 0; i < count; i++) {

        /* Get a direct pointer on the part. */
        struct part *p = &parts[pid[i]];
        struct xpart *xp = &xparts[pid[i]];

#ifdef SWIFT_DEBUG_CHECKS
        /* Is this part within the timestep? */
        if (!part_is_active(p, e)) error("Ghost applied to inactive particle");
#endif

        /* Get some useful values */
        const float h_old = p->h;
        const float h_old_dim = pow_dimension(h_old);
        const float h_old_dim_minus_one = pow_dimension_minus_one(h_old);
        float h_new;
        int has_no_neighbours = 0;

        if (p->density.wcount == 0.f) { /* No neighbours case */

          /* Flag that there were no neighbours */
          has_no_neighbours = 1;

          /* Double h and try again */
          h_new = 2.f * h_old;
        } else {

          /* Finish the density calculation */
          hydro_end_density(p, cosmo);
          chemistry_end_density(p, chemistry, cosmo);

          /* Compute one step of the Newton-Raphson scheme */
          const float n_sum = p->density.wcount * h_old_dim;
          const float n_target = hydro_eta_dim;
          const float f = n_sum - n_target;
          const float f_prime =
              p->density.wcount_dh * h_old_dim +
              hydro_dimension * p->density.wcount * h_old_dim_minus_one;

          /* Avoid floating point exception from f_prime = 0 */
          h_new = h_old - f / (f_prime + FLT_MIN);

#ifdef SWIFT_DEBUG_CHECKS
          if ((f > 0.f && h_new > h_old) || (f < 0.f && h_new < h_old))
            error(
                "Smoothing length correction not going in the right direction");
#endif

          /* Safety check: truncate to the range [ h_old/2 , 2h_old ]. */
          h_new = min(h_new, 2.f * h_old);
          h_new = max(h_new, 0.5f * h_old);
        }

        /* Check whether the particle has an inappropriate smoothing length */
        if (fabsf(h_new - h_old) > eps * h_old) {

          /* Ok, correct then */
          p->h = h_new;

          /* If below the absolute maximum, try again */
          if (p->h < hydro_h_max) {

            /* Flag for another round of fun */
            pid[redo] = pid[i];
            redo += 1;

            /* Re-initialise everything */
            hydro_init_part(p, hs);
            chemistry_init_part(p, chemistry);

            /* Off we go ! */
            continue;

          } else {

            /* Ok, this particle is a lost cause... */
            p->h = hydro_h_max;

            /* Do some damage control if no neighbours at all were found */
            if (has_no_neighbours) {
              hydro_part_has_no_neighbours(p, xp, cosmo);
              chemistry_part_has_no_neighbours(p, xp, chemistry, cosmo);
            }
          }
        }

          /* We now have a particle whose smoothing length has converged */

#ifdef EXTRA_HYDRO_LOOP

        /* As of here, particle gradient variables will be set. */
        /* The force variables are set in the extra ghost. */

        /* Compute variables required for the gradient loop */
        hydro_prepare_gradient(p, xp, cosmo);

        /* The particle gradient values are now set.  Do _NOT_
           try to read any particle density variables! */

        /* Prepare the particle for the gradient loop over neighbours */
        hydro_reset_gradient(p);

#else
        /* As of here, particle force variables will be set. */

        /* Compute variables required for the force loop */
        hydro_prepare_force(p, xp, cosmo);

        /* The particle force values are now set.  Do _NOT_
           try to read any particle density variables! */

        /* Prepare the particle for the force loop over neighbours */
        hydro_reset_acceleration(p);

#endif /* EXTRA_HYDRO_LOOP */
      }

      /* We now need to treat the particles whose smoothing length had not
       * converged again */

      /* Re-set the counter for the next loop (potentially). */
      count = redo;
      if (count > 0) {

        /* Climb up the cell hierarchy. */
        for (struct cell *finger = c; finger != NULL; finger = finger->parent) {

          /* Run through this cell's density interactions. */
          for (struct link *l = finger->density; l != NULL; l = l->next) {

#ifdef SWIFT_DEBUG_CHECKS
            if (l->t->ti_run < r->e->ti_current)
              error("Density task should have been run.");
#endif

            /* Self-interaction? */
            if (l->t->type == task_type_self)
              runner_doself_subset_branch_density(r, finger, parts, pid, count);

            /* Otherwise, pair interaction? */
            else if (l->t->type == task_type_pair) {

              /* Left or right? */
              if (l->t->ci == finger)
                runner_dopair_subset_branch_density(r, finger, parts, pid,
                                                    count, l->t->cj);
              else
                runner_dopair_subset_branch_density(r, finger, parts, pid,
                                                    count, l->t->ci);
            }

            /* Otherwise, sub-self interaction? */
            else if (l->t->type == task_type_sub_self)
              runner_dosub_subset_density(r, finger, parts, pid, count, NULL,
                                          -1, 1);

            /* Otherwise, sub-pair interaction? */
            else if (l->t->type == task_type_sub_pair) {

              /* Left or right? */
              if (l->t->ci == finger)
                runner_dosub_subset_density(r, finger, parts, pid, count,
                                            l->t->cj, -1, 1);
              else
                runner_dosub_subset_density(r, finger, parts, pid, count,
                                            l->t->ci, -1, 1);
            }
          }
        }
      }
    }

    if (count) {
      error("Smoothing length failed to converge on %i particles.", count);
    }

    /* Be clean */
    free(pid);
  }

  if (timer) TIMER_TOC(timer_do_ghost);
}

/**
 * @brief Unskip any hydro tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 */
static void runner_do_unskip_hydro(struct cell *c, struct engine *e) {

  /* Ignore empty cells. */
  if (c->count == 0) return;

  /* Skip inactive cells. */
  if (!cell_is_active_hydro(c, e)) return;

  /* Recurse */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        runner_do_unskip_hydro(cp, e);
      }
    }
  }

  /* Unskip any active tasks. */
  const int forcerebuild = cell_unskip_hydro_tasks(c, &e->sched);
  if (forcerebuild) atomic_inc(&e->forcerebuild);
}

/**
 * @brief Unskip any gravity tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 */
static void runner_do_unskip_gravity(struct cell *c, struct engine *e) {

  /* Ignore empty cells. */
  if (c->gcount == 0) return;

  /* Skip inactive cells. */
  if (!cell_is_active_gravity(c, e)) return;

  /* Recurse */
  if (c->split && c->depth < space_subdepth_grav) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        runner_do_unskip_gravity(cp, e);
      }
    }
  }

  /* Unskip any active tasks. */
  cell_unskip_gravity_tasks(c, &e->sched);
}

/**
 * @brief Mapper function to unskip active tasks.
 *
 * @param map_data An array of #cell%s.
 * @param num_elements Chunk size.
 * @param extra_data Pointer to an #engine.
 */
void runner_do_unskip_mapper(void *map_data, int num_elements,
                             void *extra_data) {

  struct engine *e = (struct engine *)extra_data;
  struct space *s = e->s;
  int *local_cells = (int *)map_data;

  for (int ind = 0; ind < num_elements; ind++) {
    struct cell *c = &s->cells_top[local_cells[ind]];
    if (c != NULL) {

      /* Hydro tasks */
      if (e->policy & engine_policy_hydro) runner_do_unskip_hydro(c, e);

      /* All gravity tasks */
      if ((e->policy & engine_policy_self_gravity) ||
          (e->policy & engine_policy_external_gravity))
        runner_do_unskip_gravity(c, e);
    }
  }
}
/**
 * @brief Drift all part in a cell.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_drift_part(struct runner *r, struct cell *c, int timer) {

  TIMER_TIC;

  cell_drift_part(c, r->e, 0);

  if (timer) TIMER_TOC(timer_drift_part);
}

/**
 * @brief Drift all gpart in a cell.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_drift_gpart(struct runner *r, struct cell *c, int timer) {

  TIMER_TIC;

  cell_drift_gpart(c, r->e, 0);

  if (timer) TIMER_TOC(timer_drift_gpart);
}

/**
 * @brief Perform the first half-kick on all the active particles in a cell.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_kick1(struct runner *r, struct cell *c, int timer) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const struct hydro_props *hydro_props = e->hydro_properties;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  struct part *restrict parts = c->parts;
  struct xpart *restrict xparts = c->xparts;
  struct gpart *restrict gparts = c->gparts;
  struct spart *restrict sparts = c->sparts;
  const int count = c->count;
  const int gcount = c->gcount;
  const int scount = c->scount;
  const integertime_t ti_current = e->ti_current;
  const double time_base = e->time_base;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_starting_hydro(c, e) && !cell_is_starting_gravity(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_kick1(r, c->progeny[k], 0);
  } else {

    /* Loop over the parts in this cell. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct part *restrict p = &parts[k];
      struct xpart *restrict xp = &xparts[k];

      /* If particle needs to be kicked */
      if (part_is_starting(p, e)) {

        const integertime_t ti_step = get_integer_timestep(p->time_bin);
        const integertime_t ti_begin =
            get_integer_time_begin(ti_current + 1, p->time_bin);

#ifdef SWIFT_DEBUG_CHECKS
        const integertime_t ti_end =
            get_integer_time_end(ti_current + 1, p->time_bin);

        if (ti_begin != ti_current)
          error(
              "Particle in wrong time-bin, ti_end=%lld, ti_begin=%lld, "
              "ti_step=%lld time_bin=%d ti_current=%lld",
              ti_end, ti_begin, ti_step, p->time_bin, ti_current);
#endif

        /* Time interval for this half-kick */
        double dt_kick_grav, dt_kick_hydro, dt_kick_therm;
        if (with_cosmology) {
          dt_kick_hydro = cosmology_get_hydro_kick_factor(
              cosmo, ti_begin, ti_begin + ti_step / 2);
          dt_kick_grav = cosmology_get_grav_kick_factor(cosmo, ti_begin,
                                                        ti_begin + ti_step / 2);
          dt_kick_therm = cosmology_get_therm_kick_factor(
              cosmo, ti_begin, ti_begin + ti_step / 2);
        } else {
          dt_kick_hydro = (ti_step / 2) * time_base;
          dt_kick_grav = (ti_step / 2) * time_base;
          dt_kick_therm = (ti_step / 2) * time_base;
        }

        /* do the kick */
        kick_part(p, xp, dt_kick_hydro, dt_kick_grav, dt_kick_therm, cosmo,
                  hydro_props, ti_begin, ti_begin + ti_step / 2);

        /* Update the accelerations to be used in the drift for hydro */
        if (p->gpart != NULL) {

          xp->a_grav[0] = p->gpart->a_grav[0];
          xp->a_grav[1] = p->gpart->a_grav[1];
          xp->a_grav[2] = p->gpart->a_grav[2];
        }
      }
    }

    /* Loop over the gparts in this cell. */
    for (int k = 0; k < gcount; k++) {

      /* Get a handle on the part. */
      struct gpart *restrict gp = &gparts[k];

      /* If the g-particle has no counterpart and needs to be kicked */
      if (gp->type == swift_type_dark_matter && gpart_is_starting(gp, e)) {

        const integertime_t ti_step = get_integer_timestep(gp->time_bin);
        const integertime_t ti_begin =
            get_integer_time_begin(ti_current + 1, gp->time_bin);

#ifdef SWIFT_DEBUG_CHECKS
        const integertime_t ti_end =
            get_integer_time_end(ti_current + 1, gp->time_bin);

        if (ti_begin != ti_current)
          error(
              "Particle in wrong time-bin, ti_end=%lld, ti_begin=%lld, "
              "ti_step=%lld time_bin=%d ti_current=%lld",
              ti_end, ti_begin, ti_step, gp->time_bin, ti_current);
#endif

        /* Time interval for this half-kick */
        double dt_kick_grav;
        if (with_cosmology) {
          dt_kick_grav = cosmology_get_grav_kick_factor(cosmo, ti_begin,
                                                        ti_begin + ti_step / 2);
        } else {
          dt_kick_grav = (ti_step / 2) * time_base;
        }

        /* do the kick */
        kick_gpart(gp, dt_kick_grav, ti_begin, ti_begin + ti_step / 2);
      }
    }

    /* Loop over the star particles in this cell. */
    for (int k = 0; k < scount; k++) {

      /* Get a handle on the s-part. */
      struct spart *restrict sp = &sparts[k];

      /* If particle needs to be kicked */
      if (spart_is_starting(sp, e)) {

        const integertime_t ti_step = get_integer_timestep(sp->time_bin);
        const integertime_t ti_begin =
            get_integer_time_begin(ti_current + 1, sp->time_bin);

#ifdef SWIFT_DEBUG_CHECKS
        const integertime_t ti_end =
            get_integer_time_end(ti_current + 1, sp->time_bin);

        if (ti_begin != ti_current)
          error(
              "Particle in wrong time-bin, ti_end=%lld, ti_begin=%lld, "
              "ti_step=%lld time_bin=%d ti_current=%lld",
              ti_end, ti_begin, ti_step, sp->time_bin, ti_current);
#endif

        /* Time interval for this half-kick */
        double dt_kick_grav;
        if (with_cosmology) {
          dt_kick_grav = cosmology_get_grav_kick_factor(cosmo, ti_begin,
                                                        ti_begin + ti_step / 2);
        } else {
          dt_kick_grav = (ti_step / 2) * time_base;
        }

        /* do the kick */
        kick_spart(sp, dt_kick_grav, ti_begin, ti_begin + ti_step / 2);
      }
    }
  }

  if (timer) TIMER_TOC(timer_kick1);
}

/**
 * @brief Perform the second half-kick on all the active particles in a cell.
 *
 * Also prepares particles to be drifted.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_kick2(struct runner *r, struct cell *c, int timer) {

  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const struct hydro_props *hydro_props = e->hydro_properties;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const int count = c->count;
  const int gcount = c->gcount;
  const int scount = c->scount;
  struct part *restrict parts = c->parts;
  struct xpart *restrict xparts = c->xparts;
  struct gpart *restrict gparts = c->gparts;
  struct spart *restrict sparts = c->sparts;
  const integertime_t ti_current = e->ti_current;
  const double time_base = e->time_base;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_hydro(c, e) && !cell_is_active_gravity(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_kick2(r, c->progeny[k], 0);
  } else {

    /* Loop over the particles in this cell. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct part *restrict p = &parts[k];
      struct xpart *restrict xp = &xparts[k];

      /* If particle needs to be kicked */
      if (part_is_active(p, e)) {

        const integertime_t ti_step = get_integer_timestep(p->time_bin);
        const integertime_t ti_begin =
            get_integer_time_begin(ti_current, p->time_bin);

#ifdef SWIFT_DEBUG_CHECKS
        if (ti_begin + ti_step != ti_current)
          error(
              "Particle in wrong time-bin, ti_begin=%lld, ti_step=%lld "
              "time_bin=%d ti_current=%lld",
              ti_begin, ti_step, p->time_bin, ti_current);
#endif
        /* Time interval for this half-kick */
        double dt_kick_grav, dt_kick_hydro, dt_kick_therm;
        if (with_cosmology) {
          dt_kick_hydro = cosmology_get_hydro_kick_factor(
              cosmo, ti_begin + ti_step / 2, ti_begin + ti_step);
          dt_kick_grav = cosmology_get_grav_kick_factor(
              cosmo, ti_begin + ti_step / 2, ti_begin + ti_step);
          dt_kick_therm = cosmology_get_therm_kick_factor(
              cosmo, ti_begin + ti_step / 2, ti_begin + ti_step);
        } else {
          dt_kick_hydro = (ti_step / 2) * time_base;
          dt_kick_grav = (ti_step / 2) * time_base;
          dt_kick_therm = (ti_step / 2) * time_base;
        }

        /* Finish the time-step with a second half-kick */
        kick_part(p, xp, dt_kick_hydro, dt_kick_grav, dt_kick_therm, cosmo,
                  hydro_props, ti_begin + ti_step / 2, ti_begin + ti_step);

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that kick and the drift are synchronized */
        if (p->ti_drift != p->ti_kick) error("Error integrating part in time.");
#endif

        /* Prepare the values to be drifted */
        hydro_reset_predicted_values(p, xp);
      }
    }

    /* Loop over the g-particles in this cell. */
    for (int k = 0; k < gcount; k++) {

      /* Get a handle on the part. */
      struct gpart *restrict gp = &gparts[k];

      /* If the g-particle has no counterpart and needs to be kicked */
      if (gp->type == swift_type_dark_matter && gpart_is_active(gp, e)) {

        const integertime_t ti_step = get_integer_timestep(gp->time_bin);
        const integertime_t ti_begin =
            get_integer_time_begin(ti_current, gp->time_bin);

#ifdef SWIFT_DEBUG_CHECKS
        if (ti_begin + ti_step != ti_current)
          error("Particle in wrong time-bin");
#endif

        /* Time interval for this half-kick */
        double dt_kick_grav;
        if (with_cosmology) {
          dt_kick_grav = cosmology_get_grav_kick_factor(
              cosmo, ti_begin + ti_step / 2, ti_begin + ti_step);
        } else {
          dt_kick_grav = (ti_step / 2) * time_base;
        }

        /* Finish the time-step with a second half-kick */
        kick_gpart(gp, dt_kick_grav, ti_begin + ti_step / 2,
                   ti_begin + ti_step);

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that kick and the drift are synchronized */
        if (gp->ti_drift != gp->ti_kick)
          error("Error integrating g-part in time.");
#endif

        /* Prepare the values to be drifted */
        gravity_reset_predicted_values(gp);
      }
    }

    /* Loop over the particles in this cell. */
    for (int k = 0; k < scount; k++) {

      /* Get a handle on the part. */
      struct spart *restrict sp = &sparts[k];

      /* If particle needs to be kicked */
      if (spart_is_active(sp, e)) {

        const integertime_t ti_step = get_integer_timestep(sp->time_bin);
        const integertime_t ti_begin =
            get_integer_time_begin(ti_current, sp->time_bin);

#ifdef SWIFT_DEBUG_CHECKS
        if (ti_begin + ti_step != ti_current)
          error("Particle in wrong time-bin");
#endif

        /* Time interval for this half-kick */
        double dt_kick_grav;
        if (with_cosmology) {
          dt_kick_grav = cosmology_get_grav_kick_factor(
              cosmo, ti_begin + ti_step / 2, ti_begin + ti_step);
        } else {
          dt_kick_grav = (ti_step / 2) * time_base;
        }

        /* Finish the time-step with a second half-kick */
        kick_spart(sp, dt_kick_grav, ti_begin + ti_step / 2,
                   ti_begin + ti_step);

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that kick and the drift are synchronized */
        if (sp->ti_drift != sp->ti_kick)
          error("Error integrating s-part in time.");
#endif

        /* Prepare the values to be drifted */
        star_reset_predicted_values(sp);
      }
    }
  }
  if (timer) TIMER_TOC(timer_kick2);
}

/**
 * @brief Computes the next time-step of all active particles in this cell
 * and update the cell's statistics.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_timestep(struct runner *r, struct cell *c, int timer) {

  const struct engine *e = r->e;
  const integertime_t ti_current = e->ti_current;
  const int count = c->count;
  const int gcount = c->gcount;
  const int scount = c->scount;
  struct part *restrict parts = c->parts;
  struct xpart *restrict xparts = c->xparts;
  struct gpart *restrict gparts = c->gparts;
  struct spart *restrict sparts = c->sparts;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_hydro(c, e) && !cell_is_active_gravity(c, e)) {
    c->updated = 0;
    c->g_updated = 0;
    c->s_updated = 0;
    return;
  }

  int updated = 0, g_updated = 0, s_updated = 0;
  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_end_max = 0,
                ti_hydro_beg_max = 0;
  integertime_t ti_gravity_end_min = max_nr_timesteps, ti_gravity_end_max = 0,
                ti_gravity_beg_max = 0;

  /* No children? */
  if (!c->split) {

    /* Loop over the particles in this cell. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct part *restrict p = &parts[k];
      struct xpart *restrict xp = &xparts[k];

      /* If particle needs updating */
      if (part_is_active(p, e)) {

#ifdef SWIFT_DEBUG_CHECKS
        /* Current end of time-step */
        const integertime_t ti_end =
            get_integer_time_end(ti_current, p->time_bin);

        if (ti_end != ti_current)
          error("Computing time-step of rogue particle.");
#endif

        /* Get new time-step */
        const integertime_t ti_new_step = get_part_timestep(p, xp, e);

        /* Update particle */
        p->time_bin = get_time_bin(ti_new_step);
        if (p->gpart != NULL) p->gpart->time_bin = p->time_bin;

        /* Number of updated particles */
        updated++;
        if (p->gpart != NULL) g_updated++;

        /* What is the next sync-point ? */
        ti_hydro_end_min = min(ti_current + ti_new_step, ti_hydro_end_min);
        ti_hydro_end_max = max(ti_current + ti_new_step, ti_hydro_end_max);

        /* What is the next starting point for this cell ? */
        ti_hydro_beg_max = max(ti_current, ti_hydro_beg_max);

        if (p->gpart != NULL) {

          /* What is the next sync-point ? */
          ti_gravity_end_min =
              min(ti_current + ti_new_step, ti_gravity_end_min);
          ti_gravity_end_max =
              max(ti_current + ti_new_step, ti_gravity_end_max);

          /* What is the next starting point for this cell ? */
          ti_gravity_beg_max = max(ti_current, ti_gravity_beg_max);
        }
      }

      else { /* part is inactive */

        const integertime_t ti_end =
            get_integer_time_end(ti_current, p->time_bin);

        const integertime_t ti_beg =
            get_integer_time_begin(ti_current + 1, p->time_bin);

        /* What is the next sync-point ? */
        ti_hydro_end_min = min(ti_end, ti_hydro_end_min);
        ti_hydro_end_max = max(ti_end, ti_hydro_end_max);

        /* What is the next starting point for this cell ? */
        ti_hydro_beg_max = max(ti_beg, ti_hydro_beg_max);

        if (p->gpart != NULL) {

          /* What is the next sync-point ? */
          ti_gravity_end_min = min(ti_end, ti_gravity_end_min);
          ti_gravity_end_max = max(ti_end, ti_gravity_end_max);

          /* What is the next starting point for this cell ? */
          ti_gravity_beg_max = max(ti_beg, ti_gravity_beg_max);
        }
      }
    }

    /* Loop over the g-particles in this cell. */
    for (int k = 0; k < gcount; k++) {

      /* Get a handle on the part. */
      struct gpart *restrict gp = &gparts[k];

      /* If the g-particle has no counterpart */
      if (gp->type == swift_type_dark_matter) {

        /* need to be updated ? */
        if (gpart_is_active(gp, e)) {

#ifdef SWIFT_DEBUG_CHECKS
          /* Current end of time-step */
          const integertime_t ti_end =
              get_integer_time_end(ti_current, gp->time_bin);

          if (ti_end != ti_current)
            error("Computing time-step of rogue particle.");
#endif

          /* Get new time-step */
          const integertime_t ti_new_step = get_gpart_timestep(gp, e);

          /* Update particle */
          gp->time_bin = get_time_bin(ti_new_step);

          /* Number of updated g-particles */
          g_updated++;

          /* What is the next sync-point ? */
          ti_gravity_end_min =
              min(ti_current + ti_new_step, ti_gravity_end_min);
          ti_gravity_end_max =
              max(ti_current + ti_new_step, ti_gravity_end_max);

          /* What is the next starting point for this cell ? */
          ti_gravity_beg_max = max(ti_current, ti_gravity_beg_max);

        } else { /* gpart is inactive */

          const integertime_t ti_end =
              get_integer_time_end(ti_current, gp->time_bin);

          /* What is the next sync-point ? */
          ti_gravity_end_min = min(ti_end, ti_gravity_end_min);
          ti_gravity_end_max = max(ti_end, ti_gravity_end_max);

          const integertime_t ti_beg =
              get_integer_time_begin(ti_current + 1, gp->time_bin);

          /* What is the next starting point for this cell ? */
          ti_gravity_beg_max = max(ti_beg, ti_gravity_beg_max);
        }
      }
    }

    /* Loop over the star particles in this cell. */
    for (int k = 0; k < scount; k++) {

      /* Get a handle on the part. */
      struct spart *restrict sp = &sparts[k];

      /* need to be updated ? */
      if (spart_is_active(sp, e)) {

#ifdef SWIFT_DEBUG_CHECKS
        /* Current end of time-step */
        const integertime_t ti_end =
            get_integer_time_end(ti_current, sp->time_bin);

        if (ti_end != ti_current)
          error("Computing time-step of rogue particle.");
#endif
        /* Get new time-step */
        const integertime_t ti_new_step = get_spart_timestep(sp, e);

        /* Update particle */
        sp->time_bin = get_time_bin(ti_new_step);
        sp->gpart->time_bin = get_time_bin(ti_new_step);

        /* Number of updated s-particles */
        s_updated++;
        g_updated++;

        /* What is the next sync-point ? */
        ti_gravity_end_min = min(ti_current + ti_new_step, ti_gravity_end_min);
        ti_gravity_end_max = max(ti_current + ti_new_step, ti_gravity_end_max);

        /* What is the next starting point for this cell ? */
        ti_gravity_beg_max = max(ti_current, ti_gravity_beg_max);

      } else { /* star particle is inactive */

        const integertime_t ti_end =
            get_integer_time_end(ti_current, sp->time_bin);

        /* What is the next sync-point ? */
        ti_gravity_end_min = min(ti_end, ti_gravity_end_min);
        ti_gravity_end_max = max(ti_end, ti_gravity_end_max);

        const integertime_t ti_beg =
            get_integer_time_begin(ti_current + 1, sp->time_bin);

        /* What is the next starting point for this cell ? */
        ti_gravity_beg_max = max(ti_beg, ti_gravity_beg_max);
      }
    }
  } else {

    /* Loop over the progeny. */
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) {
        struct cell *restrict cp = c->progeny[k];

        /* Recurse */
        runner_do_timestep(r, cp, 0);

        /* And aggregate */
        updated += cp->updated;
        g_updated += cp->g_updated;
        s_updated += cp->s_updated;
        ti_hydro_end_min = min(cp->ti_hydro_end_min, ti_hydro_end_min);
        ti_hydro_end_max = max(cp->ti_hydro_end_max, ti_hydro_end_max);
        ti_hydro_beg_max = max(cp->ti_hydro_beg_max, ti_hydro_beg_max);
        ti_gravity_end_min = min(cp->ti_gravity_end_min, ti_gravity_end_min);
        ti_gravity_end_max = max(cp->ti_gravity_end_max, ti_gravity_end_max);
        ti_gravity_beg_max = max(cp->ti_gravity_beg_max, ti_gravity_beg_max);
      }
  }

  /* Store the values. */
  c->updated = updated;
  c->g_updated = g_updated;
  c->s_updated = s_updated;
  c->ti_hydro_end_min = ti_hydro_end_min;
  c->ti_hydro_end_max = ti_hydro_end_max;
  c->ti_hydro_beg_max = ti_hydro_beg_max;
  c->ti_gravity_end_min = ti_gravity_end_min;
  c->ti_gravity_end_max = ti_gravity_end_max;
  c->ti_gravity_beg_max = ti_gravity_beg_max;

  if (timer) TIMER_TOC(timer_timestep);
}

/**
 * @brief End the force calculation of all active particles in a cell
 * by multiplying the acccelerations by the relevant constants
 *
 * @param r The #runner thread.
 * @param c The #cell.
 * @param timer Are we timing this ?
 */
void runner_do_end_force(struct runner *r, struct cell *c, int timer) {

  const struct engine *e = r->e;
  const struct space *s = e->s;
  const struct cosmology *cosmo = e->cosmology;
  const int count = c->count;
  const int gcount = c->gcount;
  const int scount = c->scount;
  struct part *restrict parts = c->parts;
  struct gpart *restrict gparts = c->gparts;
  struct spart *restrict sparts = c->sparts;
  const int periodic = s->periodic;
  const float G_newton = e->physical_constants->const_newton_G;

  TIMER_TIC;

  /* Potential normalisation in the case of periodic gravity */
  float potential_normalisation = 0.;
  if (periodic && (e->policy & engine_policy_self_gravity)) {
    const double volume = s->dim[0] * s->dim[1] * s->dim[2];
    const double r_s = e->mesh->r_s;
    potential_normalisation = 4. * M_PI * e->total_mass * r_s * r_s / volume;
  }

  /* Anything to do here? */
  if (!cell_is_active_hydro(c, e) && !cell_is_active_gravity(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_end_force(r, c->progeny[k], 0);
  } else {

    /* Loop over the gas particles in this cell. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct part *restrict p = &parts[k];

      if (part_is_active(p, e)) {

        /* Finish the force loop */
        hydro_end_force(p, cosmo);
      }
    }

    /* Loop over the g-particles in this cell. */
    for (int k = 0; k < gcount; k++) {

      /* Get a handle on the gpart. */
      struct gpart *restrict gp = &gparts[k];

      if (gpart_is_active(gp, e)) {

        /* Finish the force calculation */
        gravity_end_force(gp, G_newton, potential_normalisation, periodic);

#ifdef SWIFT_NO_GRAVITY_BELOW_ID

        /* Get the ID of the gpart */
        long long id = 0;
        if (gp->type == swift_type_gas)
          id = e->s->parts[-gp->id_or_neg_offset].id;
        else if (gp->type == swift_type_star)
          id = e->s->sparts[-gp->id_or_neg_offset].id;
        else if (gp->type == swift_type_black_hole)
          error("Unexisting type");
        else
          id = gp->id_or_neg_offset;

        /* Cancel gravity forces of these particles */
        if (id < SWIFT_NO_GRAVITY_BELOW_ID) {

          /* Don't move ! */
          gp->a_grav[0] = 0.f;
          gp->a_grav[1] = 0.f;
          gp->a_grav[2] = 0.f;
        }
#endif

#ifdef SWIFT_DEBUG_CHECKS
        if (e->policy & engine_policy_self_gravity) {

          /* Let's add a self interaction to simplify the count */
          gp->num_interacted++;

          /* Check that this gpart has interacted with all the other
           * particles (via direct or multipoles) in the box */
          if (gp->num_interacted != e->total_nr_gparts) {

            /* Get the ID of the gpart */
            long long my_id = 0;
            if (gp->type == swift_type_gas)
              my_id = e->s->parts[-gp->id_or_neg_offset].id;
            else if (gp->type == swift_type_star)
              my_id = e->s->sparts[-gp->id_or_neg_offset].id;
            else if (gp->type == swift_type_black_hole)
              error("Unexisting type");
            else
              my_id = gp->id_or_neg_offset;

            error(
                "g-particle (id=%lld, type=%s) did not interact "
                "gravitationally with all other gparts "
                "gp->num_interacted=%lld, total_gparts=%lld (local "
                "num_gparts=%zd)",
                my_id, part_type_names[gp->type], gp->num_interacted,
                e->total_nr_gparts, e->s->nr_gparts);
          }
        }
#endif
      }
    }

    /* Loop over the star particles in this cell. */
    for (int k = 0; k < scount; k++) {

      /* Get a handle on the spart. */
      struct spart *restrict sp = &sparts[k];
      if (spart_is_active(sp, e)) {

        /* Finish the force loop */
        star_end_force(sp);
      }
    }
  }

  if (timer) TIMER_TOC(timer_endforce);
}

/**
 * @brief Construct the cell properties from the received #part.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param clear_sorts Should we clear the sort flag and hence trigger a sort ?
 * @param timer Are we timing this ?
 */
void runner_do_recv_part(struct runner *r, struct cell *c, int clear_sorts,
                         int timer) {

#ifdef WITH_MPI

  const struct part *restrict parts = c->parts;
  const size_t nr_parts = c->count;
  const integertime_t ti_current = r->e->ti_current;

  TIMER_TIC;

  integertime_t ti_hydro_end_min = max_nr_timesteps;
  integertime_t ti_hydro_end_max = 0;
  timebin_t time_bin_min = num_time_bins;
  timebin_t time_bin_max = 0;
  float h_max = 0.f;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID == engine_rank) error("Updating a local cell!");
#endif

  /* Clear this cell's sorted mask. */
  if (clear_sorts) c->sorted = 0;

  /* If this cell is a leaf, collect the particle data. */
  if (!c->split) {

    /* Collect everything... */
    for (size_t k = 0; k < nr_parts; k++) {
      if (parts[k].time_bin == time_bin_inhibited) continue;
      time_bin_min = min(time_bin_min, parts[k].time_bin);
      time_bin_max = max(time_bin_max, parts[k].time_bin);
      h_max = max(h_max, parts[k].h);
    }

    /* Convert into a time */
    ti_hydro_end_min = get_integer_time_end(ti_current, time_bin_min);
    ti_hydro_end_max = get_integer_time_end(ti_current, time_bin_max);
  }

  /* Otherwise, recurse and collect. */
  else {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL && c->progeny[k]->count > 0) {
        runner_do_recv_part(r, c->progeny[k], clear_sorts, 0);
        ti_hydro_end_min =
            min(ti_hydro_end_min, c->progeny[k]->ti_hydro_end_min);
        ti_hydro_end_max =
            max(ti_hydro_end_max, c->progeny[k]->ti_hydro_end_max);
        h_max = max(h_max, c->progeny[k]->h_max);
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_hydro_end_min < ti_current)
    error(
        "Received a cell at an incorrect time c->ti_end_min=%lld, "
        "e->ti_current=%lld.",
        ti_hydro_end_min, ti_current);
#endif

  /* ... and store. */
  // c->ti_hydro_end_min = ti_hydro_end_min;
  // c->ti_hydro_end_max = ti_hydro_end_max;
  c->ti_old_part = ti_current;
  c->h_max = h_max;

  if (timer) TIMER_TOC(timer_dorecv_part);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Construct the cell properties from the received #gpart.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_recv_gpart(struct runner *r, struct cell *c, int timer) {

#ifdef WITH_MPI

  const struct gpart *restrict gparts = c->gparts;
  const size_t nr_gparts = c->gcount;
  const integertime_t ti_current = r->e->ti_current;

  TIMER_TIC;

  integertime_t ti_gravity_end_min = max_nr_timesteps;
  integertime_t ti_gravity_end_max = 0;
  timebin_t time_bin_min = num_time_bins;
  timebin_t time_bin_max = 0;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID == engine_rank) error("Updating a local cell!");
#endif

  /* If this cell is a leaf, collect the particle data. */
  if (!c->split) {

    /* Collect everything... */
    for (size_t k = 0; k < nr_gparts; k++) {
      if (gparts[k].time_bin == time_bin_inhibited) continue;
      time_bin_min = min(time_bin_min, gparts[k].time_bin);
      time_bin_max = max(time_bin_max, gparts[k].time_bin);

#ifdef SWIFT_DEBUG_CHECKS
      if (gparts[k].ti_drift != ti_current)
        error("Received un-drifted g-particle !");
#endif
    }

    /* Convert into a time */
    ti_gravity_end_min = get_integer_time_end(ti_current, time_bin_min);
    ti_gravity_end_max = get_integer_time_end(ti_current, time_bin_max);
  }

  /* Otherwise, recurse and collect. */
  else {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL && c->progeny[k]->gcount > 0) {
        runner_do_recv_gpart(r, c->progeny[k], 0);
        ti_gravity_end_min =
            min(ti_gravity_end_min, c->progeny[k]->ti_gravity_end_min);
        ti_gravity_end_max =
            max(ti_gravity_end_max, c->progeny[k]->ti_gravity_end_max);
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_gravity_end_min < ti_current)
    error(
        "Received a cell at an incorrect time c->ti_end_min=%lld, "
        "e->ti_current=%lld.",
        ti_gravity_end_min, ti_current);
#endif

  /* ... and store. */
  // c->ti_gravity_end_min = ti_gravity_end_min;
  // c->ti_gravity_end_max = ti_gravity_end_max;
  c->ti_old_gpart = ti_current;

  if (timer) TIMER_TOC(timer_dorecv_gpart);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Construct the cell properties from the received #spart.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_recv_spart(struct runner *r, struct cell *c, int timer) {

#ifdef WITH_MPI

  const struct spart *restrict sparts = c->sparts;
  const size_t nr_sparts = c->scount;
  const integertime_t ti_current = r->e->ti_current;

  TIMER_TIC;

  integertime_t ti_gravity_end_min = max_nr_timesteps;
  integertime_t ti_gravity_end_max = 0;
  timebin_t time_bin_min = num_time_bins;
  timebin_t time_bin_max = 0;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID == engine_rank) error("Updating a local cell!");
#endif

  /* If this cell is a leaf, collect the particle data. */
  if (!c->split) {

    /* Collect everything... */
    for (size_t k = 0; k < nr_sparts; k++) {
      if (sparts[k].time_bin == time_bin_inhibited) continue;
      time_bin_min = min(time_bin_min, sparts[k].time_bin);
      time_bin_max = max(time_bin_max, sparts[k].time_bin);

#ifdef SWIFT_DEBUG_CHECKS
      if (sparts[k].ti_drift != ti_current)
        error("Received un-drifted s-particle !");
#endif
    }

    /* Convert into a time */
    ti_gravity_end_min = get_integer_time_end(ti_current, time_bin_min);
    ti_gravity_end_max = get_integer_time_end(ti_current, time_bin_max);
  }

  /* Otherwise, recurse and collect. */
  else {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL && c->progeny[k]->scount > 0) {
        runner_do_recv_spart(r, c->progeny[k], 0);
        ti_gravity_end_min =
            min(ti_gravity_end_min, c->progeny[k]->ti_gravity_end_min);
        ti_gravity_end_max =
            max(ti_gravity_end_max, c->progeny[k]->ti_gravity_end_max);
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_gravity_end_min < ti_current)
    error(
        "Received a cell at an incorrect time c->ti_end_min=%lld, "
        "e->ti_current=%lld.",
        ti_gravity_end_min, ti_current);
#endif

  /* ... and store. */
  c->ti_gravity_end_min = ti_gravity_end_min;
  c->ti_gravity_end_max = ti_gravity_end_max;
  c->ti_old_gpart = ti_current;

  if (timer) TIMER_TOC(timer_dorecv_spart);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief The #runner main thread routine.
 *
 * @param data A pointer to this thread's data.
 */
void *runner_main(void *data) {

  struct runner *r = (struct runner *)data;
  struct engine *e = r->e;
  struct scheduler *sched = &e->sched;
  unsigned int seed = r->id;
  pthread_setspecific(sched->local_seed_pointer, &seed);
  /* Main loop. */
  while (1) {

    /* Wait at the barrier. */
    engine_barrier(e);

    /* Re-set the pointer to the previous task, as there is none. */
    struct task *t = NULL;
    struct task *prev = NULL;

    /* Loop while there are tasks... */
    while (1) {

      /* If there's no old task, try to get a new one. */
      if (t == NULL) {

        /* Get the task. */
        TIMER_TIC
        t = scheduler_gettask(sched, r->qid, prev);
        TIMER_TOC(timer_gettask);

        /* Did I get anything? */
        if (t == NULL) break;
      }

      /* Get the cells. */
      struct cell *ci = t->ci;
      struct cell *cj = t->cj;

#ifdef SWIFT_DEBUG_TASKS
      /* Mark the thread we run on */
      t->rid = r->cpuid;

      /* And recover the pair direction */
      if (t->type == task_type_pair || t->type == task_type_sub_pair) {
        struct cell *ci_temp = ci;
        struct cell *cj_temp = cj;
        double shift[3];
        t->sid = space_getsid(e->s, &ci_temp, &cj_temp, shift);
      } else {
        t->sid = -1;
      }
#endif

/* Check that we haven't scheduled an inactive task */
#ifdef SWIFT_DEBUG_CHECKS
      t->ti_run = e->ti_current;
#endif

      /* Different types of tasks... */
      switch (t->type) {
        case task_type_self:
          if (t->subtype == task_subtype_density)
            runner_doself1_branch_density(r, ci);
#ifdef EXTRA_HYDRO_LOOP
          else if (t->subtype == task_subtype_gradient)
            runner_doself1_branch_gradient(r, ci);
#endif
          else if (t->subtype == task_subtype_force)
            runner_doself2_branch_force(r, ci);
          else if (t->subtype == task_subtype_grav)
            runner_doself_recursive_grav(r, ci, 1);
          else if (t->subtype == task_subtype_external_grav)
            runner_do_grav_external(r, ci, 1);
          else
            error("Unknown/invalid task subtype (%d).", t->subtype);
          break;

        case task_type_pair:
          if (t->subtype == task_subtype_density)
            runner_dopair1_branch_density(r, ci, cj);
#ifdef EXTRA_HYDRO_LOOP
          else if (t->subtype == task_subtype_gradient)
            runner_dopair1_branch_gradient(r, ci, cj);
#endif
          else if (t->subtype == task_subtype_force)
            runner_dopair2_branch_force(r, ci, cj);
          else if (t->subtype == task_subtype_grav)
            runner_dopair_recursive_grav(r, ci, cj, 1);
          else
            error("Unknown/invalid task subtype (%d).", t->subtype);
          break;

        case task_type_sub_self:
          if (t->subtype == task_subtype_density)
            runner_dosub_self1_density(r, ci, 1);
#ifdef EXTRA_HYDRO_LOOP
          else if (t->subtype == task_subtype_gradient)
            runner_dosub_self1_gradient(r, ci, 1);
#endif
          else if (t->subtype == task_subtype_force)
            runner_dosub_self2_force(r, ci, 1);
          else
            error("Unknown/invalid task subtype (%d).", t->subtype);
          break;

        case task_type_sub_pair:
          if (t->subtype == task_subtype_density)
            runner_dosub_pair1_density(r, ci, cj, t->flags, 1);
#ifdef EXTRA_HYDRO_LOOP
          else if (t->subtype == task_subtype_gradient)
            runner_dosub_pair1_gradient(r, ci, cj, t->flags, 1);
#endif
          else if (t->subtype == task_subtype_force)
            runner_dosub_pair2_force(r, ci, cj, t->flags, 1);
          else
            error("Unknown/invalid task subtype (%d).", t->subtype);
          break;

        case task_type_sort:
          /* Cleanup only if any of the indices went stale. */
          runner_do_sort(r, ci, t->flags,
                         ci->dx_max_sort_old > space_maxreldx * ci->dmin, 1);
          /* Reset the sort flags as our work here is done. */
          t->flags = 0;
          break;
        case task_type_init_grav:
          runner_do_init_grav(r, ci, 1);
          break;
        case task_type_ghost:
          runner_do_ghost(r, ci, 1);
          break;
#ifdef EXTRA_HYDRO_LOOP
        case task_type_extra_ghost:
          runner_do_extra_ghost(r, ci, 1);
          break;
#endif
        case task_type_drift_part:
          runner_do_drift_part(r, ci, 1);
          break;
        case task_type_drift_gpart:
          runner_do_drift_gpart(r, ci, 1);
          break;
        case task_type_kick1:
          runner_do_kick1(r, ci, 1);
          break;
        case task_type_kick2:
          runner_do_kick2(r, ci, 1);
          break;
        case task_type_end_force:
          runner_do_end_force(r, ci, 1);
          break;
        case task_type_timestep:
          runner_do_timestep(r, ci, 1);
          break;
#ifdef WITH_MPI
        case task_type_send:
          if (t->subtype == task_subtype_tend) {
            free(t->buff);
          }
          break;
        case task_type_recv:
          if (t->subtype == task_subtype_tend) {
            cell_unpack_end_step(ci, (struct pcell_step *)t->buff);
            free(t->buff);
          } else if (t->subtype == task_subtype_xv) {
            runner_do_recv_part(r, ci, 1, 1);
          } else if (t->subtype == task_subtype_rho) {
            runner_do_recv_part(r, ci, 0, 1);
          } else if (t->subtype == task_subtype_gradient) {
            runner_do_recv_part(r, ci, 0, 1);
          } else if (t->subtype == task_subtype_gpart) {
            runner_do_recv_gpart(r, ci, 1);
          } else if (t->subtype == task_subtype_spart) {
            runner_do_recv_spart(r, ci, 1);
          } else if (t->subtype == task_subtype_multipole) {
            cell_unpack_multipoles(ci, (struct gravity_tensors *)t->buff);
            free(t->buff);
          } else {
            error("Unknown/invalid task subtype (%d).", t->subtype);
          }
          break;
#endif
        case task_type_grav_down:
          runner_do_grav_down(r, t->ci, 1);
          break;
        case task_type_grav_mesh:
          runner_do_grav_mesh(r, t->ci, 1);
          break;
        case task_type_grav_long_range:
          runner_do_grav_long_range(r, t->ci, 1);
          break;
        case task_type_grav_mm:
          runner_dopair_grav_mm_symmetric(r, t->ci, t->cj);
          break;
        case task_type_cooling:
          runner_do_cooling(r, t->ci, 1);
          break;
        case task_type_sourceterms:
          runner_do_sourceterms(r, t->ci, 1);
          break;
        default:
          error("Unknown/invalid task type (%d).", t->type);
      }

/* Mark that we have run this task on these cells */
#ifdef SWIFT_DEBUG_CHECKS
      if (ci != NULL) {
        ci->tasks_executed[t->type]++;
        ci->subtasks_executed[t->subtype]++;
      }
      if (cj != NULL) {
        cj->tasks_executed[t->type]++;
        cj->subtasks_executed[t->subtype]++;
      }
#endif

      /* We're done with this task, see if we get a next one. */
      prev = t;
      t = scheduler_done(sched, t);

    } /* main loop. */
  }

  /* Be kind, rewind. */
  return NULL;
}
