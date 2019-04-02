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
#include "entropy_floor.h"
#include "error.h"
#include "gravity.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "kick.h"
#include "logger.h"
#include "memuse.h"
#include "minmax.h"
#include "runner_doiact_vec.h"
#include "scheduler.h"
#include "sort_part.h"
#include "space.h"
#include "space_getsid.h"
#include "star_formation.h"
#include "star_formation_iact.h"
#include "stars.h"
#include "task.h"
#include "timers.h"
#include "timestep.h"
#include "timestep_limiter.h"
#include "tracers.h"

#define TASK_LOOP_DENSITY 0
#define TASK_LOOP_GRADIENT 1
#define TASK_LOOP_FORCE 2
#define TASK_LOOP_LIMITER 3
#define TASK_LOOP_FEEDBACK 4

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

/* Import the limiter loop functions. */
#define FUNCTION limiter
#define FUNCTION_TASK_LOOP TASK_LOOP_LIMITER
#include "runner_doiact.h"
#undef FUNCTION
#undef FUNCTION_TASK_LOOP

/* Import the gravity loop functions. */
#include "runner_doiact_grav.h"

/* Import the stars density loop functions. */
#define FUNCTION density
#define FUNCTION_TASK_LOOP TASK_LOOP_DENSITY
#include "runner_doiact_stars.h"
#undef FUNCTION_TASK_LOOP
#undef FUNCTION

/* Import the stars feedback loop functions. */
#define FUNCTION feedback
#define FUNCTION_TASK_LOOP TASK_LOOP_FEEDBACK
#include "runner_doiact_stars.h"
#undef FUNCTION_TASK_LOOP
#undef FUNCTION

/**
 * @brief Intermediate task after the density to check that the smoothing
 * lengths are correct.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_stars_ghost(struct runner *r, struct cell *c, int timer) {

  struct spart *restrict sparts = c->stars.parts;
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float stars_h_max = e->hydro_properties->h_max;
  const float stars_h_min = e->hydro_properties->h_min;
  const float eps = e->stars_properties->h_tolerance;
  const float stars_eta_dim =
      pow_dimension(e->stars_properties->eta_neighbours);
  const int max_smoothing_iter = e->stars_properties->max_smoothing_iterations;
  int redo = 0, scount = 0;

  TIMER_TIC;

  /* Anything to do here? */
  if (c->stars.count == 0) return;
  if (!cell_is_active_stars(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_stars_ghost(r, c->progeny[k], 0);
  } else {

    /* Init the list of active particles that have to be updated. */
    int *sid = NULL;
    float *h_0 = NULL;
    float *left = NULL;
    float *right = NULL;
    if ((sid = (int *)malloc(sizeof(int) * c->stars.count)) == NULL)
      error("Can't allocate memory for sid.");
    if ((h_0 = (float *)malloc(sizeof(float) * c->stars.count)) == NULL)
      error("Can't allocate memory for h_0.");
    if ((left = (float *)malloc(sizeof(float) * c->stars.count)) == NULL)
      error("Can't allocate memory for left.");
    if ((right = (float *)malloc(sizeof(float) * c->stars.count)) == NULL)
      error("Can't allocate memory for right.");
    for (int k = 0; k < c->stars.count; k++)
      if (spart_is_active(&sparts[k], e)) {
        sid[scount] = k;
        h_0[scount] = sparts[k].h;
        left[scount] = 0.f;
        right[scount] = stars_h_max;
        ++scount;
      }

    /* While there are particles that need to be updated... */
    for (int num_reruns = 0; scount > 0 && num_reruns < max_smoothing_iter;
         num_reruns++) {

      /* Reset the redo-count. */
      redo = 0;

      /* Loop over the remaining active parts in this cell. */
      for (int i = 0; i < scount; i++) {

        /* Get a direct pointer on the part. */
        struct spart *sp = &sparts[sid[i]];

#ifdef SWIFT_DEBUG_CHECKS
        /* Is this part within the timestep? */
        if (!spart_is_active(sp, e))
          error("Ghost applied to inactive particle");
#endif

        /* Get some useful values */
        const float h_init = h_0[i];
        const float h_old = sp->h;
        const float h_old_dim = pow_dimension(h_old);
        const float h_old_dim_minus_one = pow_dimension_minus_one(h_old);

        float h_new;
        int has_no_neighbours = 0;

        if (sp->density.wcount == 0.f) { /* No neighbours case */

          /* Flag that there were no neighbours */
          has_no_neighbours = 1;

          /* Double h and try again */
          h_new = 2.f * h_old;

        } else {

          /* Finish the density calculation */
          stars_end_density(sp, cosmo);

          /* Compute one step of the Newton-Raphson scheme */
          const float n_sum = sp->density.wcount * h_old_dim;
          const float n_target = stars_eta_dim;
          const float f = n_sum - n_target;
          const float f_prime =
              sp->density.wcount_dh * h_old_dim +
              hydro_dimension * sp->density.wcount * h_old_dim_minus_one;

          /* Improve the bisection bounds */
          if (n_sum < n_target)
            left[i] = max(left[i], h_old);
          else if (n_sum > n_target)
            right[i] = min(right[i], h_old);

#ifdef SWIFT_DEBUG_CHECKS
          /* Check the validity of the left and right bounds */
          if (left[i] > right[i])
            error("Invalid left (%e) and right (%e)", left[i], right[i]);
#endif

          /* Skip if h is already h_max and we don't have enough neighbours */
          /* Same if we are below h_min */
          if (((sp->h >= stars_h_max) && (f < 0.f)) ||
              ((sp->h <= stars_h_min) && (f > 0.f))) {

            stars_reset_feedback(sp);

            /* Ok, we are done with this particle */
            continue;
          }

          /* Normal case: Use Newton-Raphson to get a better value of h */

          /* Avoid floating point exception from f_prime = 0 */
          h_new = h_old - f / (f_prime + FLT_MIN);

          /* Be verbose about the particles that struggle to converge */
          if (num_reruns > max_smoothing_iter - 10) {

            message(
                "Smoothing length convergence problem: iter=%d p->id=%lld "
                "h_init=%12.8e h_old=%12.8e h_new=%12.8e f=%f f_prime=%f "
                "n_sum=%12.8e n_target=%12.8e left=%12.8e right=%12.8e",
                num_reruns, sp->id, h_init, h_old, h_new, f, f_prime, n_sum,
                n_target, left[i], right[i]);
          }

#ifdef SWIFT_DEBUG_CHECKS
          if ((f > 0.f && h_new > h_old) || (f < 0.f && h_new < h_old))
            error(
                "Smoothing length correction not going in the right direction");
#endif

          /* Safety check: truncate to the range [ h_old/2 , 2h_old ]. */
          h_new = min(h_new, 2.f * h_old);
          h_new = max(h_new, 0.5f * h_old);

          /* Verify that we are actually progrssing towards the answer */
          h_new = max(h_new, left[i]);
          h_new = min(h_new, right[i]);
        }

        /* Check whether the particle has an inappropriate smoothing length */
        if (fabsf(h_new - h_old) > eps * h_old) {

          /* Ok, correct then */

          /* Case where we have been oscillating around the solution */
          if ((h_new == left[i] && h_old == right[i]) ||
              (h_old == left[i] && h_new == right[i])) {

            /* Bissect the remaining interval */
            sp->h = pow_inv_dimension(
                0.5f * (pow_dimension(left[i]) + pow_dimension(right[i])));

          } else {

            /* Normal case */
            sp->h = h_new;
          }

          /* If below the absolute maximum, try again */
          if (sp->h < stars_h_max && sp->h > stars_h_min) {

            /* Flag for another round of fun */
            sid[redo] = sid[i];
            h_0[redo] = h_0[i];
            left[redo] = left[i];
            right[redo] = right[i];
            redo += 1;

            /* Re-initialise everything */
            stars_init_spart(sp);

            /* Off we go ! */
            continue;

          } else if (sp->h <= stars_h_min) {

            /* Ok, this particle is a lost cause... */
            sp->h = stars_h_min;

          } else if (sp->h >= stars_h_max) {

            /* Ok, this particle is a lost cause... */
            sp->h = stars_h_max;

            /* Do some damage control if no neighbours at all were found */
            if (has_no_neighbours) {
              stars_spart_has_no_neighbours(sp, cosmo);
            }

          } else {
            error(
                "Fundamental problem with the smoothing length iteration "
                "logic.");
          }
        }

        /* We now have a particle whose smoothing length has converged */
        stars_reset_feedback(sp);

        /* Compute the stellar evolution  */
        stars_evolve_spart(sp, e->stars_properties, cosmo);
      }

      /* We now need to treat the particles whose smoothing length had not
       * converged again */

      /* Re-set the counter for the next loop (potentially). */
      scount = redo;
      if (scount > 0) {

        /* Climb up the cell hierarchy. */
        for (struct cell *finger = c; finger != NULL; finger = finger->parent) {

          /* Run through this cell's density interactions. */
          for (struct link *l = finger->stars.density; l != NULL; l = l->next) {

#ifdef SWIFT_DEBUG_CHECKS
            if (l->t->ti_run < r->e->ti_current)
              error("Density task should have been run.");
#endif

            /* Self-interaction? */
            if (l->t->type == task_type_self)
              runner_doself_subset_branch_stars_density(r, finger, sparts, sid,
                                                        scount);

            /* Otherwise, pair interaction? */
            else if (l->t->type == task_type_pair) {

              /* Left or right? */
              if (l->t->ci == finger)
                runner_dopair_subset_branch_stars_density(
                    r, finger, sparts, sid, scount, l->t->cj);
              else
                runner_dopair_subset_branch_stars_density(
                    r, finger, sparts, sid, scount, l->t->ci);
            }

            /* Otherwise, sub-self interaction? */
            else if (l->t->type == task_type_sub_self)
              runner_dosub_subset_stars_density(r, finger, sparts, sid, scount,
                                                NULL, -1, 1);

            /* Otherwise, sub-pair interaction? */
            else if (l->t->type == task_type_sub_pair) {

              /* Left or right? */
              if (l->t->ci == finger)
                runner_dosub_subset_stars_density(r, finger, sparts, sid,
                                                  scount, l->t->cj, -1, 1);
              else
                runner_dosub_subset_stars_density(r, finger, sparts, sid,
                                                  scount, l->t->ci, -1, 1);
            }
          }
        }
      }
    }

    if (scount) {
      error("Smoothing length failed to converge on %i particles.", scount);
    }

    /* Be clean */
    free(left);
    free(right);
    free(sid);
    free(h_0);
  }

  if (timer) TIMER_TOC(timer_dostars_ghost);
}

/**
 * @brief Calculate gravity acceleration from external potential
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void runner_do_grav_external(struct runner *r, struct cell *c, int timer) {

  struct gpart *restrict gparts = c->grav.parts;
  const int gcount = c->grav.count;
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

  struct gpart *restrict gparts = c->grav.parts;
  const int gcount = c->grav.count;
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
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct entropy_floor_properties *entropy_floor_props = e->entropy_floor;
  const double time_base = e->time_base;
  const integertime_t ti_current = e->ti_current;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  const int count = c->hydro.count;

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

        double dt_cool, dt_therm;
        if (with_cosmology) {
          const integertime_t ti_step = get_integer_timestep(p->time_bin);
          const integertime_t ti_begin =
              get_integer_time_begin(ti_current - 1, p->time_bin);

          dt_cool =
              cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
          dt_therm = cosmology_get_therm_kick_factor(e->cosmology, ti_begin,
                                                     ti_begin + ti_step);

        } else {
          dt_cool = get_timestep(p->time_bin, time_base);
          dt_therm = get_timestep(p->time_bin, time_base);
        }

        /* Let's cool ! */
        cooling_cool_part(constants, us, cosmo, hydro_props,
                          entropy_floor_props, cooling_func, p, xp, dt_cool,
                          dt_therm);
      }
    }
  }

  if (timer) TIMER_TOC(timer_do_cooling);
}

/**
 *
 */
void runner_do_star_formation(struct runner *r, struct cell *c, int timer) {

  struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const struct star_formation *sf_props = e->star_formation;
  const struct phys_const *phys_const = e->physical_constants;
  const int count = c->hydro.count;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const int with_feedback = (e->policy & engine_policy_feedback);
  const struct hydro_props *restrict hydro_props = e->hydro_properties;
  const struct unit_system *restrict us = e->internal_units;
  struct cooling_function_data *restrict cooling = e->cooling_func;
  const struct entropy_floor_properties *entropy_floor = e->entropy_floor;
  const double time_base = e->time_base;
  const integertime_t ti_current = e->ti_current;
  const int current_stars_count = c->stars.count;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_hydro(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_star_formation(r, c->progeny[k], 0);
  } else {

    /* Loop over the gas particles in this cell. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct part *restrict p = &parts[k];
      struct xpart *restrict xp = &xparts[k];

      /* Only work on active particles */
      if (part_is_active(p, e)) {

        /* Is this particle star forming? */
        if (star_formation_is_star_forming(p, xp, sf_props, phys_const, cosmo,
                                           hydro_props, us, cooling,
                                           entropy_floor)) {

          /* Time-step size for this particle */
          double dt_star;
          if (with_cosmology) {
            const integertime_t ti_step = get_integer_timestep(p->time_bin);
            const integertime_t ti_begin =
                get_integer_time_begin(ti_current - 1, p->time_bin);

            dt_star =
                cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);

          } else {
            dt_star = get_timestep(p->time_bin, time_base);
          }

          /* Compute the SF rate of the particle */
          star_formation_compute_SFR(p, xp, sf_props, phys_const, cosmo,
                                     dt_star);

          /* Are we forming a star particle from this SF rate? */
          if (star_formation_should_convert_to_star(p, xp, sf_props, e,
                                                    dt_star)) {

            /* Convert the gas particle to a star particle */
            struct spart *sp = cell_convert_part_to_spart(e, c, p, xp);

            /* Copy the properties of the gas particle to the star particle */
            star_formation_copy_properties(p, xp, sp, e, sf_props, cosmo,
                                           with_cosmology);
          }

        } else { /* Are we not star-forming? */

          /* Update the particle to flag it as not star-forming */
          star_formation_update_part_not_SFR(p, xp, e, sf_props,
                                             with_cosmology);

        } /* Not Star-forming? */
      }   /* is active? */
    }     /* Loop over particles */
  }

  /* If we formed any stars, the star sorts are now invalid. We need to
   * re-compute them. */
  if (with_feedback && (c == c->hydro.super) &&
      (current_stars_count != c->stars.count)) {
    cell_clear_stars_sort_flags(c, /*is_super=*/1);
    runner_do_stars_sort(r, c, 0x1FFF, /*cleanup=*/0, /*timer=*/0);
  }

  if (timer) TIMER_TOC(timer_do_star_formation);
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

#ifdef SWIFT_DEBUG_CHECKS
/**
 * @brief Recursively checks that the flags are consistent in a cell hierarchy.
 *
 * Debugging function. Exists in two flavours: hydro & stars.
 */
#define RUNNER_CHECK_SORTS(TYPE)                                               \
  void runner_check_sorts_##TYPE(struct cell *c, int flags) {                  \
                                                                               \
    if (flags & ~c->TYPE.sorted) error("Inconsistent sort flags (downward)!"); \
    if (c->split)                                                              \
      for (int k = 0; k < 8; k++)                                              \
        if (c->progeny[k] != NULL && c->progeny[k]->TYPE.count > 0)            \
          runner_check_sorts_##TYPE(c->progeny[k], c->TYPE.sorted);            \
  }
#else
#define RUNNER_CHECK_SORTS(TYPE)                                       \
  void runner_check_sorts_##TYPE(struct cell *c, int flags) {          \
    error("Calling debugging code without debugging flag activated."); \
  }
#endif

RUNNER_CHECK_SORTS(hydro)
RUNNER_CHECK_SORTS(stars)

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
void runner_do_hydro_sort(struct runner *r, struct cell *c, int flags,
                          int cleanup, int clock) {

  struct entry *fingers[8];
  const int count = c->hydro.count;
  const struct part *parts = c->hydro.parts;
  struct xpart *xparts = c->hydro.xparts;
  float buff[8];

  TIMER_TIC;

  /* We need to do the local sorts plus whatever was requested further up. */
  flags |= c->hydro.do_sort;
  if (cleanup) {
    c->hydro.sorted = 0;
  } else {
    flags &= ~c->hydro.sorted;
  }
  if (flags == 0 && !c->hydro.do_sub_sort) return;

  /* Check that the particles have been moved to the current time */
  if (flags && !cell_are_part_drifted(c, r->e))
    error("Sorting un-drifted cell c->nodeID=%d", c->nodeID);

#ifdef SWIFT_DEBUG_CHECKS
  /* Make sure the sort flags are consistent (downward). */
  runner_check_sorts_hydro(c, c->hydro.sorted);

  /* Make sure the sort flags are consistent (upard). */
  for (struct cell *finger = c->parent; finger != NULL;
       finger = finger->parent) {
    if (finger->hydro.sorted & ~c->hydro.sorted)
      error("Inconsistent sort flags (upward).");
  }

  /* Update the sort timer which represents the last time the sorts
     were re-set. */
  if (c->hydro.sorted == 0) c->hydro.ti_sort = r->e->ti_current;
#endif

  /* Allocate memory for sorting. */
  cell_malloc_hydro_sorts(c, flags);

  /* Does this cell have any progeny? */
  if (c->split) {

    /* Fill in the gaps within the progeny. */
    float dx_max_sort = 0.0f;
    float dx_max_sort_old = 0.0f;
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL && c->progeny[k]->hydro.count > 0) {
        /* Only propagate cleanup if the progeny is stale. */
        runner_do_hydro_sort(r, c->progeny[k], flags,
                             cleanup && (c->progeny[k]->hydro.dx_max_sort_old >
                                         space_maxreldx * c->progeny[k]->dmin),
                             0);
        dx_max_sort = max(dx_max_sort, c->progeny[k]->hydro.dx_max_sort);
        dx_max_sort_old =
            max(dx_max_sort_old, c->progeny[k]->hydro.dx_max_sort_old);
      }
    }
    c->hydro.dx_max_sort = dx_max_sort;
    c->hydro.dx_max_sort_old = dx_max_sort_old;

    /* Loop over the 13 different sort arrays. */
    for (int j = 0; j < 13; j++) {

      /* Has this sort array been flagged? */
      if (!(flags & (1 << j))) continue;

      /* Init the particle index offsets. */
      int off[8];
      off[0] = 0;
      for (int k = 1; k < 8; k++)
        if (c->progeny[k - 1] != NULL)
          off[k] = off[k - 1] + c->progeny[k - 1]->hydro.count;
        else
          off[k] = off[k - 1];

      /* Init the entries and indices. */
      int inds[8];
      for (int k = 0; k < 8; k++) {
        inds[k] = k;
        if (c->progeny[k] != NULL && c->progeny[k]->hydro.count > 0) {
          fingers[k] = c->progeny[k]->hydro.sort[j];
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
      struct entry *finger = c->hydro.sort[j];
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
      c->hydro.sort[j][count].d = FLT_MAX;
      c->hydro.sort[j][count].i = 0;

      /* Mark as sorted. */
      atomic_or(&c->hydro.sorted, 1 << j);

    } /* loop over sort arrays. */

  } /* progeny? */

  /* Otherwise, just sort. */
  else {

    /* Reset the sort distance */
    if (c->hydro.sorted == 0) {
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
      c->hydro.dx_max_sort_old = 0.f;
      c->hydro.dx_max_sort = 0.f;
    }

    /* Fill the sort array. */
    for (int k = 0; k < count; k++) {
      const double px[3] = {parts[k].x[0], parts[k].x[1], parts[k].x[2]};
      for (int j = 0; j < 13; j++)
        if (flags & (1 << j)) {
          c->hydro.sort[j][k].i = k;
          c->hydro.sort[j][k].d = px[0] * runner_shift[j][0] +
                                  px[1] * runner_shift[j][1] +
                                  px[2] * runner_shift[j][2];
        }
    }

    /* Add the sentinel and sort. */
    for (int j = 0; j < 13; j++)
      if (flags & (1 << j)) {
        c->hydro.sort[j][count].d = FLT_MAX;
        c->hydro.sort[j][count].i = 0;
        runner_do_sort_ascending(c->hydro.sort[j], count);
        atomic_or(&c->hydro.sorted, 1 << j);
      }
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify the sorting. */
  for (int j = 0; j < 13; j++) {
    if (!(flags & (1 << j))) continue;
    struct entry *finger = c->hydro.sort[j];
    for (int k = 1; k < count; k++) {
      if (finger[k].d < finger[k - 1].d)
        error("Sorting failed, ascending array.");
      if (finger[k].i >= count) error("Sorting failed, indices borked.");
    }
  }

  /* Make sure the sort flags are consistent (downward). */
  runner_check_sorts_hydro(c, flags);

  /* Make sure the sort flags are consistent (upward). */
  for (struct cell *finger = c->parent; finger != NULL;
       finger = finger->parent) {
    if (finger->hydro.sorted & ~c->hydro.sorted)
      error("Inconsistent sort flags.");
  }
#endif

  /* Clear the cell's sort flags. */
  c->hydro.do_sort = 0;
  c->hydro.do_sub_sort = 0;
  c->hydro.requires_sorts = 0;

  if (clock) TIMER_TOC(timer_dosort);
}

/**
 * @brief Sort the stars particles in the given cell along all cardinal
 * directions.
 *
 * @param r The #runner.
 * @param c The #cell.
 * @param flags Cell flag.
 * @param cleanup If true, re-build the sorts for the selected flags instead
 *        of just adding them.
 * @param clock Flag indicating whether to record the timing or not, needed
 *      for recursive calls.
 */
void runner_do_stars_sort(struct runner *r, struct cell *c, int flags,
                          int cleanup, int clock) {

  struct entry *fingers[8];
  const int count = c->stars.count;
  struct spart *sparts = c->stars.parts;
  float buff[8];

  TIMER_TIC;

  /* We need to do the local sorts plus whatever was requested further up. */
  flags |= c->stars.do_sort;
  if (cleanup) {
    c->stars.sorted = 0;
  } else {
    flags &= ~c->stars.sorted;
  }
  if (flags == 0 && !c->stars.do_sub_sort) return;

  /* Check that the particles have been moved to the current time */
  if (flags && !cell_are_spart_drifted(c, r->e)) {
    error("Sorting un-drifted cell c->nodeID=%d", c->nodeID);
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Make sure the sort flags are consistent (downward). */
  runner_check_sorts_stars(c, c->stars.sorted);

  /* Make sure the sort flags are consistent (upward). */
  for (struct cell *finger = c->parent; finger != NULL;
       finger = finger->parent) {
    if (finger->stars.sorted & ~c->stars.sorted)
      error("Inconsistent sort flags (upward).");
  }

  /* Update the sort timer which represents the last time the sorts
     were re-set. */
  if (c->stars.sorted == 0) c->stars.ti_sort = r->e->ti_current;
#endif

  /* start by allocating the entry arrays in the requested dimensions. */
  cell_malloc_stars_sorts(c, flags);

  /* Does this cell have any progeny? */
  if (c->split) {

    /* Fill in the gaps within the progeny. */
    float dx_max_sort = 0.0f;
    float dx_max_sort_old = 0.0f;
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL && c->progeny[k]->stars.count > 0) {
        /* Only propagate cleanup if the progeny is stale. */
        const int cleanup_prog =
            cleanup && (c->progeny[k]->stars.dx_max_sort_old >
                        space_maxreldx * c->progeny[k]->dmin);
        runner_do_stars_sort(r, c->progeny[k], flags, cleanup_prog, 0);
        dx_max_sort = max(dx_max_sort, c->progeny[k]->stars.dx_max_sort);
        dx_max_sort_old =
            max(dx_max_sort_old, c->progeny[k]->stars.dx_max_sort_old);
      }
    }
    c->stars.dx_max_sort = dx_max_sort;
    c->stars.dx_max_sort_old = dx_max_sort_old;

    /* Loop over the 13 different sort arrays. */
    for (int j = 0; j < 13; j++) {

      /* Has this sort array been flagged? */
      if (!(flags & (1 << j))) continue;

      /* Init the particle index offsets. */
      int off[8];
      off[0] = 0;
      for (int k = 1; k < 8; k++)
        if (c->progeny[k - 1] != NULL)
          off[k] = off[k - 1] + c->progeny[k - 1]->stars.count;
        else
          off[k] = off[k - 1];

      /* Init the entries and indices. */
      int inds[8];
      for (int k = 0; k < 8; k++) {
        inds[k] = k;
        if (c->progeny[k] != NULL && c->progeny[k]->stars.count > 0) {
          fingers[k] = c->progeny[k]->stars.sort[j];
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
      struct entry *finger = c->stars.sort[j];
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
      c->stars.sort[j][count].d = FLT_MAX;
      c->stars.sort[j][count].i = 0;

      /* Mark as sorted. */
      atomic_or(&c->stars.sorted, 1 << j);

    } /* loop over sort arrays. */

  } /* progeny? */

  /* Otherwise, just sort. */
  else {

    /* Reset the sort distance */
    if (c->stars.sorted == 0) {

      /* And the individual sort distances if we are a local cell */
      for (int k = 0; k < count; k++) {
        sparts[k].x_diff_sort[0] = 0.0f;
        sparts[k].x_diff_sort[1] = 0.0f;
        sparts[k].x_diff_sort[2] = 0.0f;
      }
      c->stars.dx_max_sort_old = 0.f;
      c->stars.dx_max_sort = 0.f;
    }

    /* Fill the sort array. */
    for (int k = 0; k < count; k++) {
      const double px[3] = {sparts[k].x[0], sparts[k].x[1], sparts[k].x[2]};
      for (int j = 0; j < 13; j++)
        if (flags & (1 << j)) {
          c->stars.sort[j][k].i = k;
          c->stars.sort[j][k].d = px[0] * runner_shift[j][0] +
                                  px[1] * runner_shift[j][1] +
                                  px[2] * runner_shift[j][2];
        }
    }

    /* Add the sentinel and sort. */
    for (int j = 0; j < 13; j++)
      if (flags & (1 << j)) {
        c->stars.sort[j][count].d = FLT_MAX;
        c->stars.sort[j][count].i = 0;
        runner_do_sort_ascending(c->stars.sort[j], count);
        atomic_or(&c->stars.sorted, 1 << j);
      }
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* Verify the sorting. */
  for (int j = 0; j < 13; j++) {
    if (!(flags & (1 << j))) continue;
    struct entry *finger = c->stars.sort[j];
    for (int k = 1; k < count; k++) {
      if (finger[k].d < finger[k - 1].d)
        error("Sorting failed, ascending array.");
      if (finger[k].i >= count) error("Sorting failed, indices borked.");
    }
  }

  /* Make sure the sort flags are consistent (downward). */
  runner_check_sorts_stars(c, flags);

  /* Make sure the sort flags are consistent (upward). */
  for (struct cell *finger = c->parent; finger != NULL;
       finger = finger->parent) {
    if (finger->stars.sorted & ~c->stars.sorted)
      error("Inconsistent sort flags.");
  }
#endif

  /* Clear the cell's sort flags. */
  c->stars.do_sort = 0;
  c->stars.do_sub_sort = 0;
  c->stars.requires_sorts = 0;

  if (clock) TIMER_TOC(timer_do_stars_sort);
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
  gravity_field_tensors_init(&c->grav.multipole->pot, e->ti_current);

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

  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  const int count = c->hydro.count;
  const struct engine *e = r->e;
  const integertime_t ti_end = e->ti_current;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const double time_base = e->time_base;
  const struct cosmology *cosmo = e->cosmology;
  const struct hydro_props *hydro_props = e->hydro_properties;

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

        /* Calculate the time-step for passing to hydro_prepare_force.
         * This is the physical time between the start and end of the time-step
         * without any scale-factor powers. */
        double dt_alpha;
        if (with_cosmology) {
          const integertime_t ti_step = get_integer_timestep(p->time_bin);
          dt_alpha = cosmology_get_delta_time(cosmo, ti_end - ti_step, ti_end);
        } else {
          dt_alpha = get_timestep(p->time_bin, time_base);
        }

        /* Compute variables required for the force loop */
        hydro_prepare_force(p, xp, cosmo, hydro_props, dt_alpha);

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

  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  const struct engine *e = r->e;
  const struct space *s = e->s;
  const struct hydro_space *hs = &s->hs;
  const struct cosmology *cosmo = e->cosmology;
  const struct chemistry_global_data *chemistry = e->chemistry;
  const struct star_formation *star_formation = e->star_formation;
  const float hydro_h_max = e->hydro_properties->h_max;
  const float hydro_h_min = e->hydro_properties->h_min;
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

    /* Init the list of active particles that have to be updated and their
     * current smoothing lengths. */
    int *pid = NULL;
    float *h_0 = NULL;
    float *left = NULL;
    float *right = NULL;
    if ((pid = (int *)malloc(sizeof(int) * c->hydro.count)) == NULL)
      error("Can't allocate memory for pid.");
    if ((h_0 = (float *)malloc(sizeof(float) * c->hydro.count)) == NULL)
      error("Can't allocate memory for h_0.");
    if ((left = (float *)malloc(sizeof(float) * c->hydro.count)) == NULL)
      error("Can't allocate memory for left.");
    if ((right = (float *)malloc(sizeof(float) * c->hydro.count)) == NULL)
      error("Can't allocate memory for right.");
    for (int k = 0; k < c->hydro.count; k++)
      if (part_is_active(&parts[k], e)) {
        pid[count] = k;
        h_0[count] = parts[k].h;
        left[count] = 0.f;
        right[count] = hydro_h_max;
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
        const float h_init = h_0[i];
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
          star_formation_end_density(p, star_formation, cosmo);

          /* Compute one step of the Newton-Raphson scheme */
          const float n_sum = p->density.wcount * h_old_dim;
          const float n_target = hydro_eta_dim;
          const float f = n_sum - n_target;
          const float f_prime =
              p->density.wcount_dh * h_old_dim +
              hydro_dimension * p->density.wcount * h_old_dim_minus_one;

          /* Improve the bisection bounds */
          if (n_sum < n_target)
            left[i] = max(left[i], h_old);
          else if (n_sum > n_target)
            right[i] = min(right[i], h_old);

#ifdef SWIFT_DEBUG_CHECKS
          /* Check the validity of the left and right bounds */
          if (left[i] > right[i])
            error("Invalid left (%e) and right (%e)", left[i], right[i]);
#endif

          /* Skip if h is already h_max and we don't have enough neighbours */
          /* Same if we are below h_min */
          if (((p->h >= hydro_h_max) && (f < 0.f)) ||
              ((p->h <= hydro_h_min) && (f > 0.f))) {

          /* We have a particle whose smoothing length is already set (wants
           * to be larger but has already hit the maximum OR wants to be
           * smaller but has already reached the minimum). So, just tidy up as
           * if the smoothing length had converged correctly  */

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
            const struct hydro_props *hydro_props = e->hydro_properties;

            /* Calculate the time-step for passing to hydro_prepare_force, used
             * for the evolution of alpha factors (i.e. those involved in the
             * artificial viscosity and thermal conduction terms) */
            const int with_cosmology = (e->policy & engine_policy_cosmology);
            const double time_base = e->time_base;
            const integertime_t ti_end = e->ti_current;
            double dt_alpha;

            if (with_cosmology) {
              const integertime_t ti_step = get_integer_timestep(p->time_bin);
              dt_alpha =
                  cosmology_get_delta_time(cosmo, ti_end - ti_step, ti_end);
            } else {
              dt_alpha = get_timestep(p->time_bin, time_base);
            }

            /* As of here, particle force variables will be set. */

            /* Compute variables required for the force loop */
            hydro_prepare_force(p, xp, cosmo, hydro_props, dt_alpha);

            /* The particle force values are now set.  Do _NOT_
               try to read any particle density variables! */

            /* Prepare the particle for the force loop over neighbours */
            hydro_reset_acceleration(p);

#endif /* EXTRA_HYDRO_LOOP */

            /* Ok, we are done with this particle */
            continue;
          }

          /* Normal case: Use Newton-Raphson to get a better value of h */

          /* Avoid floating point exception from f_prime = 0 */
          h_new = h_old - f / (f_prime + FLT_MIN);

          /* Be verbose about the particles that struggle to converge */
          if (num_reruns > max_smoothing_iter - 10) {

            message(
                "Smoothing length convergence problem: iter=%d p->id=%lld "
                "h_init=%12.8e h_old=%12.8e h_new=%12.8e f=%f f_prime=%f "
                "n_sum=%12.8e n_target=%12.8e left=%12.8e right=%12.8e",
                num_reruns, p->id, h_init, h_old, h_new, f, f_prime, n_sum,
                n_target, left[i], right[i]);
          }

#ifdef SWIFT_DEBUG_CHECKS
          if ((f > 0.f && h_new > h_old) || (f < 0.f && h_new < h_old))
            error(
                "Smoothing length correction not going in the right direction");
#endif

          /* Safety check: truncate to the range [ h_old/2 , 2h_old ]. */
          h_new = min(h_new, 2.f * h_old);
          h_new = max(h_new, 0.5f * h_old);

          /* Verify that we are actually progrssing towards the answer */
          h_new = max(h_new, left[i]);
          h_new = min(h_new, right[i]);
        }

        /* Check whether the particle has an inappropriate smoothing length */
        if (fabsf(h_new - h_old) > eps * h_old) {

          /* Ok, correct then */

          /* Case where we have been oscillating around the solution */
          if ((h_new == left[i] && h_old == right[i]) ||
              (h_old == left[i] && h_new == right[i])) {

            /* Bissect the remaining interval */
            p->h = pow_inv_dimension(
                0.5f * (pow_dimension(left[i]) + pow_dimension(right[i])));

          } else {

            /* Normal case */
            p->h = h_new;
          }

          /* If within the allowed range, try again */
          if (p->h < hydro_h_max && p->h > hydro_h_min) {

            /* Flag for another round of fun */
            pid[redo] = pid[i];
            h_0[redo] = h_0[i];
            left[redo] = left[i];
            right[redo] = right[i];
            redo += 1;

            /* Re-initialise everything */
            hydro_init_part(p, hs);
            chemistry_init_part(p, chemistry);
            star_formation_init_part(p, star_formation);

            /* Off we go ! */
            continue;

          } else if (p->h <= hydro_h_min) {

            /* Ok, this particle is a lost cause... */
            p->h = hydro_h_min;

          } else if (p->h >= hydro_h_max) {

            /* Ok, this particle is a lost cause... */
            p->h = hydro_h_max;

            /* Do some damage control if no neighbours at all were found */
            if (has_no_neighbours) {
              hydro_part_has_no_neighbours(p, xp, cosmo);
              chemistry_part_has_no_neighbours(p, xp, chemistry, cosmo);
              star_formation_part_has_no_neighbours(p, xp, star_formation,
                                                    cosmo);
            }

          } else {
            error(
                "Fundamental problem with the smoothing length iteration "
                "logic.");
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
        const struct hydro_props *hydro_props = e->hydro_properties;

        /* Calculate the time-step for passing to hydro_prepare_force, used for
         * the evolution of alpha factors (i.e. those involved in the artificial
         * viscosity and thermal conduction terms) */
        const int with_cosmology = (e->policy & engine_policy_cosmology);
        const integertime_t ti_end = e->ti_current;
        const double time_base = e->time_base;
        double dt_alpha;

        if (with_cosmology) {
          const integertime_t ti_step = get_integer_timestep(p->time_bin);
          dt_alpha = cosmology_get_delta_time(cosmo, ti_end - ti_step, ti_end);
        } else {
          dt_alpha = get_timestep(p->time_bin, time_base);
        }

        /* As of here, particle force variables will be set. */

        /* Compute variables required for the force loop */
        hydro_prepare_force(p, xp, cosmo, hydro_props, dt_alpha);

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
          for (struct link *l = finger->hydro.density; l != NULL; l = l->next) {

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
    free(left);
    free(right);
    free(pid);
    free(h_0);
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
  if (c->hydro.count == 0) return;

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
 * @brief Unskip any stars tasks associated with active cells.
 *
 * @param c The cell.
 * @param e The engine.
 */
static void runner_do_unskip_stars(struct cell *c, struct engine *e) {

  /* Ignore empty cells. */
  if (c->stars.count == 0) return;

  /* Skip inactive cells. */
  if (!cell_is_active_stars(c, e)) return;

  /* Recurse */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];
        runner_do_unskip_stars(cp, e);
      }
    }
  }

  /* Unskip any active tasks. */
  const int forcerebuild = cell_unskip_stars_tasks(c, &e->sched);
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
  if (c->grav.count == 0) return;

  /* Skip inactive cells. */
  if (!cell_is_active_gravity(c, e)) return;

  /* Recurse */
  if (c->split && ((c->maxdepth - c->depth) >= space_subdepth_diff_grav)) {
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

      /* Stars tasks */
      if (e->policy & engine_policy_stars) runner_do_unskip_stars(c, e);
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
 * @brief Drift all spart in a cell.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_drift_spart(struct runner *r, struct cell *c, int timer) {

  TIMER_TIC;

  cell_drift_spart(c, r->e, 0);

  if (timer) TIMER_TOC(timer_drift_spart);
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
  const struct entropy_floor_properties *entropy_floor = e->entropy_floor;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  struct gpart *restrict gparts = c->grav.parts;
  struct spart *restrict sparts = c->stars.parts;
  const int count = c->hydro.count;
  const int gcount = c->grav.count;
  const int scount = c->stars.count;
  const integertime_t ti_current = e->ti_current;
  const double time_base = e->time_base;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_starting_hydro(c, e) && !cell_is_starting_gravity(c, e) &&
      !cell_is_starting_stars(c, e))
    return;

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

#ifdef SWIFT_DEBUG_CHECKS
        if (p->wakeup == time_bin_awake)
          error("Woken-up particle that has not been processed in kick1");
#endif

        /* Skip particles that have been woken up and treated by the limiter. */
        if (p->wakeup != time_bin_not_awake) continue;

        const integertime_t ti_step = get_integer_timestep(p->time_bin);
        const integertime_t ti_begin =
            get_integer_time_begin(ti_current + 1, p->time_bin);

#ifdef SWIFT_DEBUG_CHECKS
        const integertime_t ti_end = ti_begin + ti_step;

        if (ti_begin != ti_current)
          error(
              "Particle in wrong time-bin, ti_end=%lld, ti_begin=%lld, "
              "ti_step=%lld time_bin=%d wakeup=%d ti_current=%lld",
              ti_end, ti_begin, ti_step, p->time_bin, p->wakeup, ti_current);
#endif

        /* Time interval for this half-kick */
        double dt_kick_grav, dt_kick_hydro, dt_kick_therm, dt_kick_corr;
        if (with_cosmology) {
          dt_kick_hydro = cosmology_get_hydro_kick_factor(
              cosmo, ti_begin, ti_begin + ti_step / 2);
          dt_kick_grav = cosmology_get_grav_kick_factor(cosmo, ti_begin,
                                                        ti_begin + ti_step / 2);
          dt_kick_therm = cosmology_get_therm_kick_factor(
              cosmo, ti_begin, ti_begin + ti_step / 2);
          dt_kick_corr = cosmology_get_corr_kick_factor(cosmo, ti_begin,
                                                        ti_begin + ti_step / 2);
        } else {
          dt_kick_hydro = (ti_step / 2) * time_base;
          dt_kick_grav = (ti_step / 2) * time_base;
          dt_kick_therm = (ti_step / 2) * time_base;
          dt_kick_corr = (ti_step / 2) * time_base;
        }

        /* do the kick */
        kick_part(p, xp, dt_kick_hydro, dt_kick_grav, dt_kick_therm,
                  dt_kick_corr, cosmo, hydro_props, entropy_floor, ti_begin,
                  ti_begin + ti_step / 2);

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

    /* Loop over the stars particles in this cell. */
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
  const struct entropy_floor_properties *entropy_floor = e->entropy_floor;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const int count = c->hydro.count;
  const int gcount = c->grav.count;
  const int scount = c->stars.count;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  struct gpart *restrict gparts = c->grav.parts;
  struct spart *restrict sparts = c->stars.parts;
  const integertime_t ti_current = e->ti_current;
  const double time_base = e->time_base;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_hydro(c, e) && !cell_is_active_gravity(c, e) &&
      !cell_is_active_stars(c, e))
    return;

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

        integertime_t ti_begin, ti_end, ti_step;

#ifdef SWIFT_DEBUG_CHECKS
        if (p->wakeup == time_bin_awake)
          error("Woken-up particle that has not been processed in kick1");
#endif

        if (p->wakeup == time_bin_not_awake) {

          /* Time-step from a regular kick */
          ti_step = get_integer_timestep(p->time_bin);
          ti_begin = get_integer_time_begin(ti_current, p->time_bin);
          ti_end = ti_begin + ti_step;

        } else {

          /* Time-step that follows a wake-up call */
          ti_begin = get_integer_time_begin(ti_current, p->wakeup);
          ti_end = get_integer_time_end(ti_current, p->time_bin);
          ti_step = ti_end - ti_begin;

          /* Reset the flag. Everything is back to normal from now on. */
          p->wakeup = time_bin_awake;
        }

#ifdef SWIFT_DEBUG_CHECKS
        if (ti_begin + ti_step != ti_current)
          error(
              "Particle in wrong time-bin, ti_begin=%lld, ti_step=%lld "
              "time_bin=%d wakeup=%d ti_current=%lld",
              ti_begin, ti_step, p->time_bin, p->wakeup, ti_current);
#endif
        /* Time interval for this half-kick */
        double dt_kick_grav, dt_kick_hydro, dt_kick_therm, dt_kick_corr;
        if (with_cosmology) {
          dt_kick_hydro = cosmology_get_hydro_kick_factor(
              cosmo, ti_begin + ti_step / 2, ti_end);
          dt_kick_grav = cosmology_get_grav_kick_factor(
              cosmo, ti_begin + ti_step / 2, ti_end);
          dt_kick_therm = cosmology_get_therm_kick_factor(
              cosmo, ti_begin + ti_step / 2, ti_end);
          dt_kick_corr = cosmology_get_corr_kick_factor(
              cosmo, ti_begin + ti_step / 2, ti_end);
        } else {
          dt_kick_hydro = (ti_end - (ti_begin + ti_step / 2)) * time_base;
          dt_kick_grav = (ti_end - (ti_begin + ti_step / 2)) * time_base;
          dt_kick_therm = (ti_end - (ti_begin + ti_step / 2)) * time_base;
          dt_kick_corr = (ti_end - (ti_begin + ti_step / 2)) * time_base;
        }

        /* Finish the time-step with a second half-kick */
        kick_part(p, xp, dt_kick_hydro, dt_kick_grav, dt_kick_therm,
                  dt_kick_corr, cosmo, hydro_props, entropy_floor,
                  ti_begin + ti_step / 2, ti_end);

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
        stars_reset_predicted_values(sp);
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
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const int count = c->hydro.count;
  const int gcount = c->grav.count;
  const int scount = c->stars.count;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  struct gpart *restrict gparts = c->grav.parts;
  struct spart *restrict sparts = c->stars.parts;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_hydro(c, e) && !cell_is_active_gravity(c, e) &&
      !cell_is_active_stars(c, e)) {
    c->hydro.updated = 0;
    c->grav.updated = 0;
    c->stars.updated = 0;
    return;
  }

  int updated = 0, g_updated = 0, s_updated = 0;
  int inhibited = 0, g_inhibited = 0, s_inhibited = 0;
  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_end_max = 0,
                ti_hydro_beg_max = 0;
  integertime_t ti_gravity_end_min = max_nr_timesteps, ti_gravity_end_max = 0,
                ti_gravity_beg_max = 0;
  integertime_t ti_stars_end_min = max_nr_timesteps, ti_stars_end_max = 0,
                ti_stars_beg_max = 0;

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

        /* Update the tracers properties */
        tracers_after_timestep(p, xp, e->internal_units, e->physical_constants,
                               with_cosmology, e->cosmology,
                               e->hydro_properties, e->cooling_func, e->time);

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

        /* Count the number of inhibited particles */
        if (part_is_inhibited(p, e)) inhibited++;

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

          /* Count the number of inhibited particles */
          if (gpart_is_inhibited(gp, e)) g_inhibited++;

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

        ti_stars_end_min = min(ti_current + ti_new_step, ti_stars_end_min);
        ti_stars_end_max = max(ti_current + ti_new_step, ti_stars_end_max);
        ti_gravity_end_min = min(ti_current + ti_new_step, ti_gravity_end_min);
        ti_gravity_end_max = max(ti_current + ti_new_step, ti_gravity_end_max);

        /* What is the next starting point for this cell ? */
        ti_stars_beg_max = max(ti_current, ti_stars_beg_max);
        ti_gravity_beg_max = max(ti_current, ti_gravity_beg_max);

        /* star particle is inactive but not inhibited */
      } else {

        /* Count the number of inhibited particles */
        if (spart_is_inhibited(sp, e)) ++s_inhibited;

        const integertime_t ti_end =
            get_integer_time_end(ti_current, sp->time_bin);

        const integertime_t ti_beg =
            get_integer_time_begin(ti_current + 1, sp->time_bin);

        ti_stars_end_min = min(ti_end, ti_stars_end_min);
        ti_stars_end_max = max(ti_end, ti_stars_end_max);
        ti_gravity_end_min = min(ti_end, ti_gravity_end_min);
        ti_gravity_end_max = max(ti_end, ti_gravity_end_max);

        /* What is the next starting point for this cell ? */
        ti_stars_beg_max = max(ti_beg, ti_stars_beg_max);
        ti_gravity_beg_max = max(ti_beg, ti_gravity_beg_max);
      }
    }
  } else {

    /* Loop over the progeny. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *restrict cp = c->progeny[k];

        /* Recurse */
        runner_do_timestep(r, cp, 0);

        /* And aggregate */
        updated += cp->hydro.updated;
        g_updated += cp->grav.updated;
        s_updated += cp->stars.updated;
        inhibited += cp->hydro.inhibited;
        g_inhibited += cp->grav.inhibited;
        s_inhibited += cp->stars.inhibited;
        ti_hydro_end_min = min(cp->hydro.ti_end_min, ti_hydro_end_min);
        ti_hydro_end_max = max(cp->hydro.ti_end_max, ti_hydro_end_max);
        ti_hydro_beg_max = max(cp->hydro.ti_beg_max, ti_hydro_beg_max);
        ti_gravity_end_min = min(cp->grav.ti_end_min, ti_gravity_end_min);
        ti_gravity_end_max = max(cp->grav.ti_end_max, ti_gravity_end_max);
        ti_gravity_beg_max = max(cp->grav.ti_beg_max, ti_gravity_beg_max);
        ti_stars_end_min = min(cp->stars.ti_end_min, ti_stars_end_min);
        ti_stars_end_max = max(cp->grav.ti_end_max, ti_stars_end_max);
        ti_stars_beg_max = max(cp->grav.ti_beg_max, ti_stars_beg_max);
      }
    }
  }

  /* Store the values. */
  c->hydro.updated = updated;
  c->grav.updated = g_updated;
  c->stars.updated = s_updated;
  c->hydro.inhibited = inhibited;
  c->grav.inhibited = g_inhibited;
  c->stars.inhibited = s_inhibited;
  c->hydro.ti_end_min = ti_hydro_end_min;
  c->hydro.ti_end_max = ti_hydro_end_max;
  c->hydro.ti_beg_max = ti_hydro_beg_max;
  c->grav.ti_end_min = ti_gravity_end_min;
  c->grav.ti_end_max = ti_gravity_end_max;
  c->grav.ti_beg_max = ti_gravity_beg_max;
  c->stars.ti_end_min = ti_stars_end_min;
  c->stars.ti_end_max = ti_stars_end_max;
  c->stars.ti_beg_max = ti_stars_beg_max;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->hydro.ti_end_min == e->ti_current &&
      c->hydro.ti_end_min < max_nr_timesteps)
    error("End of next hydro step is current time!");
  if (c->grav.ti_end_min == e->ti_current &&
      c->grav.ti_end_min < max_nr_timesteps)
    error("End of next gravity step is current time!");
  if (c->stars.ti_end_min == e->ti_current &&
      c->stars.ti_end_min < max_nr_timesteps)
    error("End of next stars step is current time!");
#endif

  if (timer) TIMER_TOC(timer_timestep);
}

/**
 * @brief Apply the time-step limiter to all awaken particles in a cell
 * hierarchy.
 *
 * @param r The task #runner.
 * @param c The #cell.
 * @param force Limit the particles irrespective of the #cell flags.
 * @param timer Are we timing this ?
 */
void runner_do_limiter(struct runner *r, struct cell *c, int force, int timer) {

  const struct engine *e = r->e;
  const integertime_t ti_current = e->ti_current;
  const int count = c->hydro.count;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we only limit local cells. */
  if (c->nodeID != engine_rank) error("Limiting dt of a foreign cell is nope.");
#endif

  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_end_max = 0,
                ti_hydro_beg_max = 0;
  integertime_t ti_gravity_end_min = max_nr_timesteps, ti_gravity_end_max = 0,
                ti_gravity_beg_max = 0;

  /* Limit irrespective of cell flags? */
  force |= c->hydro.do_limiter;

  /* Early abort? */
  if (c->hydro.count == 0) {

    /* Clear the limiter flags. */
    c->hydro.do_limiter = 0;
    c->hydro.do_sub_limiter = 0;
    return;
  }

  /* Loop over the progeny ? */
  if (c->split && (force || c->hydro.do_sub_limiter)) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *restrict cp = c->progeny[k];

        /* Recurse */
        runner_do_limiter(r, cp, force, 0);

        /* And aggregate */
        ti_hydro_end_min = min(cp->hydro.ti_end_min, ti_hydro_end_min);
        ti_hydro_end_max = max(cp->hydro.ti_end_max, ti_hydro_end_max);
        ti_hydro_beg_max = max(cp->hydro.ti_beg_max, ti_hydro_beg_max);
        ti_gravity_end_min = min(cp->grav.ti_end_min, ti_gravity_end_min);
        ti_gravity_end_max = max(cp->grav.ti_end_max, ti_gravity_end_max);
        ti_gravity_beg_max = max(cp->grav.ti_beg_max, ti_gravity_beg_max);
      }
    }

    /* Store the updated values */
    c->hydro.ti_end_min = min(c->hydro.ti_end_min, ti_hydro_end_min);
    c->hydro.ti_end_max = max(c->hydro.ti_end_max, ti_hydro_end_max);
    c->hydro.ti_beg_max = max(c->hydro.ti_beg_max, ti_hydro_beg_max);
    c->grav.ti_end_min = min(c->grav.ti_end_min, ti_gravity_end_min);
    c->grav.ti_end_max = max(c->grav.ti_end_max, ti_gravity_end_max);
    c->grav.ti_beg_max = max(c->grav.ti_beg_max, ti_gravity_beg_max);

  } else if (!c->split && force) {

    ti_hydro_end_min = c->hydro.ti_end_min;
    ti_hydro_end_max = c->hydro.ti_end_max;
    ti_hydro_beg_max = c->hydro.ti_beg_max;
    ti_gravity_end_min = c->grav.ti_end_min;
    ti_gravity_end_max = c->grav.ti_end_max;
    ti_gravity_beg_max = c->grav.ti_beg_max;

    /* Loop over the gas particles in this cell. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct part *restrict p = &parts[k];
      struct xpart *restrict xp = &xparts[k];

      /* Avoid inhibited particles */
      if (part_is_inhibited(p, e)) continue;

      /* If the particle will be active no need to wake it up */
      if (part_is_active(p, e) && p->wakeup != time_bin_not_awake)
        p->wakeup = time_bin_not_awake;

      /* Bip, bip, bip... wake-up time */
      if (p->wakeup == time_bin_awake) {

        /* Apply the limiter and get the new time-step size */
        const integertime_t ti_new_step = timestep_limit_part(p, xp, e);

        /* What is the next sync-point ? */
        ti_hydro_end_min = min(ti_current + ti_new_step, ti_hydro_end_min);
        ti_hydro_end_max = max(ti_current + ti_new_step, ti_hydro_end_max);

        /* What is the next starting point for this cell ? */
        ti_hydro_beg_max = max(ti_current, ti_hydro_beg_max);

        /* Also limit the gpart counter-part */
        if (p->gpart != NULL) {

          /* Register the time-bin */
          p->gpart->time_bin = p->time_bin;

          /* What is the next sync-point ? */
          ti_gravity_end_min =
              min(ti_current + ti_new_step, ti_gravity_end_min);
          ti_gravity_end_max =
              max(ti_current + ti_new_step, ti_gravity_end_max);

          /* What is the next starting point for this cell ? */
          ti_gravity_beg_max = max(ti_current, ti_gravity_beg_max);
        }
      }
    }

    /* Store the updated values */
    c->hydro.ti_end_min = min(c->hydro.ti_end_min, ti_hydro_end_min);
    c->hydro.ti_end_max = max(c->hydro.ti_end_max, ti_hydro_end_max);
    c->hydro.ti_beg_max = max(c->hydro.ti_beg_max, ti_hydro_beg_max);
    c->grav.ti_end_min = min(c->grav.ti_end_min, ti_gravity_end_min);
    c->grav.ti_end_max = max(c->grav.ti_end_max, ti_gravity_end_max);
    c->grav.ti_beg_max = max(c->grav.ti_beg_max, ti_gravity_beg_max);
  }

  /* Clear the limiter flags. */
  c->hydro.do_limiter = 0;
  c->hydro.do_sub_limiter = 0;

  if (timer) TIMER_TOC(timer_do_limiter);
}

/**
 * @brief End the hydro force calculation of all active particles in a cell
 * by multiplying the acccelerations by the relevant constants
 *
 * @param r The #runner thread.
 * @param c The #cell.
 * @param timer Are we timing this ?
 */
void runner_do_end_hydro_force(struct runner *r, struct cell *c, int timer) {

  const struct engine *e = r->e;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_hydro(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_end_hydro_force(r, c->progeny[k], 0);
  } else {

    const struct cosmology *cosmo = e->cosmology;
    const int count = c->hydro.count;
    struct part *restrict parts = c->hydro.parts;

    /* Loop over the gas particles in this cell. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct part *restrict p = &parts[k];

      if (part_is_active(p, e)) {

        /* Finish the force loop */
        hydro_end_force(p, cosmo);
      }
    }
  }

  if (timer) TIMER_TOC(timer_end_hydro_force);
}

/**
 * @brief End the gravity force calculation of all active particles in a cell
 * by multiplying the acccelerations by the relevant constants
 *
 * @param r The #runner thread.
 * @param c The #cell.
 * @param timer Are we timing this ?
 */
void runner_do_end_grav_force(struct runner *r, struct cell *c, int timer) {

  const struct engine *e = r->e;

  TIMER_TIC;

  /* Anything to do here? */
  if (!cell_is_active_gravity(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_end_grav_force(r, c->progeny[k], 0);
  } else {

    const struct space *s = e->s;
    const int periodic = s->periodic;
    const float G_newton = e->physical_constants->const_newton_G;

    /* Potential normalisation in the case of periodic gravity */
    float potential_normalisation = 0.;
    if (periodic && (e->policy & engine_policy_self_gravity)) {
      const double volume = s->dim[0] * s->dim[1] * s->dim[2];
      const double r_s = e->mesh->r_s;
      potential_normalisation = 4. * M_PI * e->total_mass * r_s * r_s / volume;
    }

    const int gcount = c->grav.count;
    struct gpart *restrict gparts = c->grav.parts;

    /* Loop over the g-particles in this cell. */
    for (int k = 0; k < gcount; k++) {

      /* Get a handle on the gpart. */
      struct gpart *restrict gp = &gparts[k];

      if (gpart_is_active(gp, e)) {

        /* Finish the force calculation */
        gravity_end_force(gp, G_newton, potential_normalisation, periodic);

#ifdef SWIFT_MAKE_GRAVITY_GLASS

        /* Negate the gravity forces */
        gp->a_grav[0] *= -1.f;
        gp->a_grav[1] *= -1.f;
        gp->a_grav[2] *= -1.f;
#endif

#ifdef SWIFT_NO_GRAVITY_BELOW_ID

        /* Get the ID of the gpart */
        long long id = 0;
        if (gp->type == swift_type_gas)
          id = e->s->parts[-gp->id_or_neg_offset].id;
        else if (gp->type == swift_type_stars)
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
          if (gp->num_interacted !=
              e->total_nr_gparts - e->count_inhibited_gparts) {

            /* Get the ID of the gpart */
            long long my_id = 0;
            if (gp->type == swift_type_gas)
              my_id = e->s->parts[-gp->id_or_neg_offset].id;
            else if (gp->type == swift_type_stars)
              my_id = e->s->sparts[-gp->id_or_neg_offset].id;
            else if (gp->type == swift_type_black_hole)
              error("Unexisting type");
            else
              my_id = gp->id_or_neg_offset;

            error(
                "g-particle (id=%lld, type=%s) did not interact "
                "gravitationally with all other gparts "
                "gp->num_interacted=%lld, total_gparts=%lld (local "
                "num_gparts=%zd inhibited_gparts=%lld)",
                my_id, part_type_names[gp->type], gp->num_interacted,
                e->total_nr_gparts, e->s->nr_gparts, e->count_inhibited_gparts);
          }
        }
#endif
      }
    }
  }
  if (timer) TIMER_TOC(timer_end_grav_force);
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

  const struct part *restrict parts = c->hydro.parts;
  const size_t nr_parts = c->hydro.count;
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
  if (clear_sorts) c->hydro.sorted = 0;

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
      if (c->progeny[k] != NULL && c->progeny[k]->hydro.count > 0) {
        runner_do_recv_part(r, c->progeny[k], clear_sorts, 0);
        ti_hydro_end_min =
            min(ti_hydro_end_min, c->progeny[k]->hydro.ti_end_min);
        ti_hydro_end_max =
            max(ti_hydro_end_max, c->progeny[k]->hydro.ti_end_max);
        h_max = max(h_max, c->progeny[k]->hydro.h_max);
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
  // c->hydro.ti_end_min = ti_hydro_end_min;
  // c->hydro.ti_end_max = ti_hydro_end_max;
  c->hydro.ti_old_part = ti_current;
  c->hydro.h_max = h_max;

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

  const struct gpart *restrict gparts = c->grav.parts;
  const size_t nr_gparts = c->grav.count;
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
    }

    /* Convert into a time */
    ti_gravity_end_min = get_integer_time_end(ti_current, time_bin_min);
    ti_gravity_end_max = get_integer_time_end(ti_current, time_bin_max);
  }

  /* Otherwise, recurse and collect. */
  else {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL && c->progeny[k]->grav.count > 0) {
        runner_do_recv_gpart(r, c->progeny[k], 0);
        ti_gravity_end_min =
            min(ti_gravity_end_min, c->progeny[k]->grav.ti_end_min);
        ti_gravity_end_max =
            max(ti_gravity_end_max, c->progeny[k]->grav.ti_end_max);
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
  // c->grav.ti_end_min = ti_gravity_end_min;
  // c->grav.ti_end_max = ti_gravity_end_max;
  c->grav.ti_old_part = ti_current;

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
 * @param clear_sorts Should we clear the sort flag and hence trigger a sort ?
 * @param timer Are we timing this ?
 */
void runner_do_recv_spart(struct runner *r, struct cell *c, int clear_sorts,
                          int timer) {

#ifdef WITH_MPI

  struct spart *restrict sparts = c->stars.parts;
  const size_t nr_sparts = c->stars.count;
  const integertime_t ti_current = r->e->ti_current;

  TIMER_TIC;

  integertime_t ti_stars_end_min = max_nr_timesteps;
  integertime_t ti_stars_end_max = 0;
  timebin_t time_bin_min = num_time_bins;
  timebin_t time_bin_max = 0;
  float h_max = 0.f;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID == engine_rank) error("Updating a local cell!");
#endif

  /* Clear this cell's sorted mask. */
  if (clear_sorts) c->stars.sorted = 0;

  /* If this cell is a leaf, collect the particle data. */
  if (!c->split) {

    /* Collect everything... */
    for (size_t k = 0; k < nr_sparts; k++) {
#ifdef DEBUG_INTERACTIONS_STARS
      sparts[k].num_ngb_force = 0;
#endif
      if (sparts[k].time_bin == time_bin_inhibited) continue;
      time_bin_min = min(time_bin_min, sparts[k].time_bin);
      time_bin_max = max(time_bin_max, sparts[k].time_bin);
      h_max = max(h_max, sparts[k].h);
    }

    /* Convert into a time */
    ti_stars_end_min = get_integer_time_end(ti_current, time_bin_min);
    ti_stars_end_max = get_integer_time_end(ti_current, time_bin_max);
  }

  /* Otherwise, recurse and collect. */
  else {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL && c->progeny[k]->stars.count > 0) {
        runner_do_recv_spart(r, c->progeny[k], clear_sorts, 0);
        ti_stars_end_min =
            min(ti_stars_end_min, c->progeny[k]->stars.ti_end_min);
        ti_stars_end_max =
            max(ti_stars_end_max, c->progeny[k]->grav.ti_end_max);
        h_max = max(h_max, c->progeny[k]->stars.h_max);
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_stars_end_min < ti_current)
    error(
        "Received a cell at an incorrect time c->ti_end_min=%lld, "
        "e->ti_current=%lld.",
        ti_stars_end_min, ti_current);
#endif

  /* ... and store. */
  // c->grav.ti_end_min = ti_gravity_end_min;
  // c->grav.ti_end_max = ti_gravity_end_max;
  c->stars.ti_old_part = ti_current;
  c->stars.h_max = h_max;

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
          else if (t->subtype == task_subtype_limiter)
            runner_doself2_branch_limiter(r, ci);
          else if (t->subtype == task_subtype_grav)
            runner_doself_recursive_grav(r, ci, 1);
          else if (t->subtype == task_subtype_external_grav)
            runner_do_grav_external(r, ci, 1);
          else if (t->subtype == task_subtype_stars_density)
            runner_doself_branch_stars_density(r, ci);
          else if (t->subtype == task_subtype_stars_feedback)
            runner_doself_branch_stars_feedback(r, ci);
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
          else if (t->subtype == task_subtype_limiter)
            runner_dopair2_branch_limiter(r, ci, cj);
          else if (t->subtype == task_subtype_grav)
            runner_dopair_recursive_grav(r, ci, cj, 1);
          else if (t->subtype == task_subtype_stars_density)
            runner_dopair_branch_stars_density(r, ci, cj);
          else if (t->subtype == task_subtype_stars_feedback)
            runner_dopair_branch_stars_feedback(r, ci, cj);
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
          else if (t->subtype == task_subtype_limiter)
            runner_dosub_self2_limiter(r, ci, 1);
          else if (t->subtype == task_subtype_stars_density)
            runner_dosub_self_stars_density(r, ci, 1);
          else if (t->subtype == task_subtype_stars_feedback)
            runner_dosub_self_stars_feedback(r, ci, 1);
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
          else if (t->subtype == task_subtype_limiter)
            runner_dosub_pair2_limiter(r, ci, cj, t->flags, 1);
          else if (t->subtype == task_subtype_stars_density)
            runner_dosub_pair_stars_density(r, ci, cj, t->flags, 1);
          else if (t->subtype == task_subtype_stars_feedback)
            runner_dosub_pair_stars_feedback(r, ci, cj, t->flags, 1);
          else
            error("Unknown/invalid task subtype (%d).", t->subtype);
          break;

        case task_type_sort:
          /* Cleanup only if any of the indices went stale. */
          runner_do_hydro_sort(
              r, ci, t->flags,
              ci->hydro.dx_max_sort_old > space_maxreldx * ci->dmin, 1);
          /* Reset the sort flags as our work here is done. */
          t->flags = 0;
          break;
        case task_type_stars_sort:
          /* Cleanup only if any of the indices went stale. */
          runner_do_stars_sort(
              r, ci, t->flags,
              ci->stars.dx_max_sort_old > space_maxreldx * ci->dmin, 1);
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
        case task_type_stars_ghost:
          runner_do_stars_ghost(r, ci, 1);
          break;
        case task_type_drift_part:
          runner_do_drift_part(r, ci, 1);
          break;
        case task_type_drift_spart:
          runner_do_drift_spart(r, ci, 1);
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
        case task_type_end_hydro_force:
          runner_do_end_hydro_force(r, ci, 1);
          break;
        case task_type_end_grav_force:
          runner_do_end_grav_force(r, ci, 1);
          break;
        case task_type_logger:
          runner_do_logger(r, ci, 1);
          break;
        case task_type_timestep:
          runner_do_timestep(r, ci, 1);
          break;
        case task_type_timestep_limiter:
          runner_do_limiter(r, ci, 0, 1);
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
          } else if (t->subtype == task_subtype_limiter) {
            runner_do_recv_part(r, ci, 0, 1);
          } else if (t->subtype == task_subtype_gpart) {
            runner_do_recv_gpart(r, ci, 1);
          } else if (t->subtype == task_subtype_spart) {
            runner_do_recv_spart(r, ci, 1, 1);
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
          runner_dopair_grav_mm_progenies(r, t->flags, t->ci, t->cj);
          break;
        case task_type_cooling:
          runner_do_cooling(r, t->ci, 1);
          break;
        case task_type_star_formation:
          runner_do_star_formation(r, t->ci, 1);
          break;
        case task_type_fof_self:
          runner_do_fof_self(r, t->ci, 1);
          break;
        case task_type_fof_pair:
          runner_do_fof_pair(r, t->ci, t->cj, 1);
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

/**
 * @brief Write the required particles through the logger.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_logger(struct runner *r, struct cell *c, int timer) {

#ifdef WITH_LOGGER
  TIMER_TIC;

  const struct engine *e = r->e;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  const int count = c->hydro.count;

  /* Anything to do here? */
  if (!cell_is_starting_hydro(c, e) && !cell_is_starting_gravity(c, e)) return;

  /* Recurse? Avoid spending too much time in useless cells. */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) runner_do_logger(r, c->progeny[k], 0);
  } else {

    /* Loop over the parts in this cell. */
    for (int k = 0; k < count; k++) {

      /* Get a handle on the part. */
      struct part *restrict p = &parts[k];
      struct xpart *restrict xp = &xparts[k];

      /* If particle needs to be log */
      /* This is the same function than part_is_active, except for
       * debugging checks */
      if (part_is_starting(p, e)) {

        if (logger_should_write(&xp->logger_data, e->logger)) {
          /* Write particle */
          /* Currently writing everything, should adapt it through time */
          logger_log_part(e->logger, p,
                          logger_mask_data[logger_x].mask |
                              logger_mask_data[logger_v].mask |
                              logger_mask_data[logger_a].mask |
                              logger_mask_data[logger_u].mask |
                              logger_mask_data[logger_h].mask |
                              logger_mask_data[logger_rho].mask |
                              logger_mask_data[logger_consts].mask,
                          &xp->logger_data.last_offset);

          /* Set counter back to zero */
          xp->logger_data.steps_since_last_output = 0;
        } else
          /* Update counter */
          xp->logger_data.steps_since_last_output += 1;
      }
    }
  }

  if (c->grav.count > 0) error("gparts not implemented");

  if (c->stars.count > 0) error("sparts not implemented");

  if (timer) TIMER_TOC(timer_logger);

#else
  error("Logger disabled, please enable it during configuration");
#endif
}

/**
 * @brief Recursively search for FOF groups in a single cell.
 *
 * @param r runner task
 * @param c cell
 * @param timer 1 if the time is to be recorded.
 */
void runner_do_fof_self(struct runner *r, struct cell *c, int timer) {

  TIMER_TIC;

  struct space *s = r->e->s;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double search_r2 = s->fof_data.l_x2;

  rec_fof_search_self(c, s, dim, search_r2);

  if (timer) TIMER_TOC(timer_fof_self);
}

/**
 * @brief Recursively search for FOF groups between a pair of cells.
 *
 * @param r runner task
 * @param ci cell i
 * @param cj cell j
 * @param timer 1 if the time is to be recorded.
 */
void runner_do_fof_pair(struct runner *r, struct cell *ci, struct cell *cj,
                        int timer) {

  TIMER_TIC;

  struct space *s = r->e->s;
  const double dim[3] = {s->dim[0], s->dim[1], s->dim[2]};
  const double search_r2 = s->fof_data.l_x2;

  rec_fof_search_pair(ci, cj, s, dim, search_r2);

  if (timer) TIMER_TOC(timer_fof_pair);
}
