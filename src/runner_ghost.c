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
#include "runner.h"

/* Local headers. */
#include "active.h"
#include "black_holes.h"
#include "cell.h"
#include "dark_matter.h"
#include "engine.h"
#include "feedback.h"
#include "pressure_floor.h"
#include "pressure_floor_iact.h"
#include "space_getsid.h"
#include "star_formation.h"
#include "stars.h"
#include "timers.h"
#include "timestep_limiter.h"
#include "tracers.h"

/* Import the dark matter functions. */
#include "runner_doiact_dark_matter.h"

/* Import the density loop functions. */
#define FUNCTION density
#define FUNCTION_TASK_LOOP TASK_LOOP_DENSITY
#include "runner_doiact_hydro.h"
#undef FUNCTION
#undef FUNCTION_TASK_LOOP

/* Import the stars density loop functions. */
#define FUNCTION density
#define FUNCTION_TASK_LOOP TASK_LOOP_DENSITY
#include "runner_doiact_stars.h"
#undef FUNCTION_TASK_LOOP
#undef FUNCTION

/* Import the black hole density loop functions. */
#define FUNCTION density
#define FUNCTION_TASK_LOOP TASK_LOOP_DENSITY
#include "runner_doiact_black_holes.h"
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
  const struct unit_system *us = e->internal_units;
  const struct phys_const *phys_const = e->physical_constants;
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const struct cosmology *cosmo = e->cosmology;
  const struct feedback_props *feedback_props = e->feedback_props;
  const float stars_h_max = e->hydro_properties->h_max;
  const float stars_h_min = e->hydro_properties->h_min;
  const float eps = e->stars_properties->h_tolerance;
  const float stars_eta_dim =
      pow_dimension(e->stars_properties->eta_neighbours);
  const int max_smoothing_iter = e->stars_properties->max_smoothing_iterations;
  int redo = 0, scount = 0;

  /* Running value of the maximal smoothing length */
  double h_max = c->stars.h_max;

  TIMER_TIC;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID != e->nodeID)
    error("Running the star ghost on a foreign node!");
#endif

  /* Anything to do here? */
  if (c->stars.count == 0) return;
  if (!cell_is_active_stars(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        runner_do_stars_ghost(r, c->progeny[k], 0);

        /* Update h_max */
        h_max = max(h_max, c->progeny[k]->stars.h_max);
      }
    }
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
      if (spart_is_active(&sparts[k], e) &&
          feedback_is_active(&sparts[k], e->time, cosmo, with_cosmology)) {
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

          /* Skip if h is already h_max and we don't have enough neighbours
           */
          /* Same if we are below h_min */
          if (((sp->h >= stars_h_max) && (f < 0.f)) ||
              ((sp->h <= stars_h_min) && (f > 0.f))) {

            stars_reset_feedback(sp);

            /* Only do feedback if stars have a reasonable birth time */
            if (feedback_is_active(sp, e->time, cosmo, with_cosmology)) {

              const integertime_t ti_step = get_integer_timestep(sp->time_bin);
              const integertime_t ti_begin =
                  get_integer_time_begin(e->ti_current - 1, sp->time_bin);

              /* Get particle time-step */
              double dt_star;
              if (with_cosmology) {
                dt_star = cosmology_get_delta_time(e->cosmology, ti_begin,
                                                   ti_begin + ti_step);
              } else {
                dt_star = get_timestep(sp->time_bin, e->time_base);
              }

              /* Calculate age of the star at current time */
              double star_age_end_of_step;
              if (with_cosmology) {
                star_age_end_of_step =
                    cosmology_get_delta_time_from_scale_factors(
                        cosmo, (double)sp->birth_scale_factor, cosmo->a);
              } else {
                star_age_end_of_step = e->time - (double)sp->birth_time;
              }

              /* Has this star been around for a while ? */
              if (star_age_end_of_step > 0.) {

                /* Get the length of the enrichment time-step */
                const double dt_enrichment = feedback_get_enrichment_timestep(
                    sp, with_cosmology, cosmo, e->time, dt_star);
                const double star_age_beg_of_step =
                    star_age_end_of_step - dt_enrichment;

                /* Compute the stellar evolution  */
                feedback_evolve_spart(sp, feedback_props, cosmo, us, phys_const,
                                      star_age_beg_of_step, dt_enrichment,
                                      e->time, ti_begin, with_cosmology);
              } else {

                /* Reset the feedback fields of the star particle */
                feedback_reset_feedback(sp, feedback_props);
              }
            } else {

              feedback_reset_feedback(sp, feedback_props);
            }

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

          /* Safety check: truncate to the range [ h_old/2 , 2h_old ]. */
          h_new = min(h_new, 2.f * h_old);
          h_new = max(h_new, 0.5f * h_old);

          /* Verify that we are actually progrssing towards the answer */
          h_new = max(h_new, left[i]);
          h_new = min(h_new, right[i]);
        }

        /* Check whether the particle has an inappropriate smoothing length
         */
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
            feedback_init_spart(sp);

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

        /* Check if h_max has increased */
        h_max = max(h_max, sp->h);

        stars_reset_feedback(sp);

        /* Only do feedback if stars have a reasonable birth time */
        if (feedback_is_active(sp, e->time, cosmo, with_cosmology)) {

          const integertime_t ti_step = get_integer_timestep(sp->time_bin);
          const integertime_t ti_begin =
              get_integer_time_begin(e->ti_current - 1, sp->time_bin);

          /* Get particle time-step */
          double dt_star;
          if (with_cosmology) {
            dt_star = cosmology_get_delta_time(e->cosmology, ti_begin,
                                               ti_begin + ti_step);
          } else {
            dt_star = get_timestep(sp->time_bin, e->time_base);
          }

          /* Calculate age of the star at current time */
          double star_age_end_of_step;
          if (with_cosmology) {
            star_age_end_of_step = cosmology_get_delta_time_from_scale_factors(
                cosmo, (double)sp->birth_scale_factor, cosmo->a);
          } else {
            star_age_end_of_step = e->time - (double)sp->birth_time;
          }

          /* Has this star been around for a while ? */
          if (star_age_end_of_step > 0.) {

            /* Get the length of the enrichment time-step */
            const double dt_enrichment = feedback_get_enrichment_timestep(
                sp, with_cosmology, cosmo, e->time, dt_star);
            const double star_age_beg_of_step =
                star_age_end_of_step - dt_enrichment;

            /* Compute the stellar evolution  */
            feedback_evolve_spart(sp, feedback_props, cosmo, us, phys_const,
                                  star_age_beg_of_step, dt_enrichment, e->time,
                                  ti_begin, with_cosmology);
          } else {

            /* Reset the feedback fields of the star particle */
            feedback_reset_feedback(sp, feedback_props);
          }
        } else {

          /* Reset the feedback fields of the star particle */
          feedback_reset_feedback(sp, feedback_props);
        }
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
                                                NULL, 1);

            /* Otherwise, sub-pair interaction? */
            else if (l->t->type == task_type_sub_pair) {

              /* Left or right? */
              if (l->t->ci == finger)
                runner_dosub_subset_stars_density(r, finger, sparts, sid,
                                                  scount, l->t->cj, 1);
              else
                runner_dosub_subset_stars_density(r, finger, sparts, sid,
                                                  scount, l->t->ci, 1);
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

  /* Update h_max */
  c->stars.h_max = h_max;

  /* The ghost may not always be at the top level.
   * Therefore we need to update h_max between the super- and top-levels */
  if (c->stars.ghost) {
    for (struct cell *tmp = c->parent; tmp != NULL; tmp = tmp->parent) {
      atomic_max_f(&tmp->stars.h_max, h_max);
    }
  }

  if (timer) TIMER_TOC(timer_do_stars_ghost);
}

/**
 * @brief Intermediate task after the density to check that the smoothing
 * lengths are correct.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_black_holes_density_ghost(struct runner *r, struct cell *c,
                                         int timer) {

  struct bpart *restrict bparts = c->black_holes.parts;
  const struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  const float black_holes_h_max = e->hydro_properties->h_max;
  const float black_holes_h_min = e->hydro_properties->h_min;
  const float eps = e->black_holes_properties->h_tolerance;
  const float black_holes_eta_dim =
      pow_dimension(e->black_holes_properties->eta_neighbours);
  const int max_smoothing_iter = e->hydro_properties->max_smoothing_iterations;
  int redo = 0, bcount = 0;

  /* Running value of the maximal smoothing length */
  double h_max = c->black_holes.h_max;

  TIMER_TIC;

  /* Anything to do here? */
  if (c->black_holes.count == 0) return;
  if (!cell_is_active_black_holes(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        runner_do_black_holes_density_ghost(r, c->progeny[k], 0);

        /* Update h_max */
        h_max = max(h_max, c->progeny[k]->black_holes.h_max);
      }
    }
  } else {

    /* Init the list of active particles that have to be updated. */
    int *sid = NULL;
    float *h_0 = NULL;
    float *left = NULL;
    float *right = NULL;
    if ((sid = (int *)malloc(sizeof(int) * c->black_holes.count)) == NULL)
      error("Can't allocate memory for sid.");
    if ((h_0 = (float *)malloc(sizeof(float) * c->black_holes.count)) == NULL)
      error("Can't allocate memory for h_0.");
    if ((left = (float *)malloc(sizeof(float) * c->black_holes.count)) == NULL)
      error("Can't allocate memory for left.");
    if ((right = (float *)malloc(sizeof(float) * c->black_holes.count)) == NULL)
      error("Can't allocate memory for right.");
    for (int k = 0; k < c->black_holes.count; k++)
      if (bpart_is_active(&bparts[k], e)) {
        sid[bcount] = k;
        h_0[bcount] = bparts[k].h;
        left[bcount] = 0.f;
        right[bcount] = black_holes_h_max;
        ++bcount;
      }

    /* While there are particles that need to be updated... */
    for (int num_reruns = 0; bcount > 0 && num_reruns < max_smoothing_iter;
         num_reruns++) {

      /* Reset the redo-count. */
      redo = 0;

      /* Loop over the remaining active parts in this cell. */
      for (int i = 0; i < bcount; i++) {

        /* Get a direct pointer on the part. */
        struct bpart *bp = &bparts[sid[i]];

#ifdef SWIFT_DEBUG_CHECKS
        /* Is this part within the timestep? */
        if (!bpart_is_active(bp, e))
          error("Ghost applied to inactive particle");
#endif

        /* Get some useful values */
        const float h_init = h_0[i];
        const float h_old = bp->h;
        const float h_old_dim = pow_dimension(h_old);
        const float h_old_dim_minus_one = pow_dimension_minus_one(h_old);

        float h_new;
        int has_no_neighbours = 0;

        if (bp->density.wcount == 0.f) { /* No neighbours case */

          /* Flag that there were no neighbours */
          has_no_neighbours = 1;

          /* Double h and try again */
          h_new = 2.f * h_old;

        } else {

          /* Finish the density calculation */
          black_holes_end_density(bp, cosmo);

          /* Compute one step of the Newton-Raphson scheme */
          const float n_sum = bp->density.wcount * h_old_dim;
          const float n_target = black_holes_eta_dim;
          const float f = n_sum - n_target;
          const float f_prime =
              bp->density.wcount_dh * h_old_dim +
              hydro_dimension * bp->density.wcount * h_old_dim_minus_one;

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

          /* Skip if h is already h_max and we don't have enough neighbours
           */
          /* Same if we are below h_min */
          if (((bp->h >= black_holes_h_max) && (f < 0.f)) ||
              ((bp->h <= black_holes_h_min) && (f > 0.f))) {

            black_holes_reset_feedback(bp);

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
                num_reruns, bp->id, h_init, h_old, h_new, f, f_prime, n_sum,
                n_target, left[i], right[i]);
          }

          /* Safety check: truncate to the range [ h_old/2 , 2h_old ]. */
          h_new = min(h_new, 2.f * h_old);
          h_new = max(h_new, 0.5f * h_old);

          /* Verify that we are actually progrssing towards the answer */
          h_new = max(h_new, left[i]);
          h_new = min(h_new, right[i]);
        }

        /* Check whether the particle has an inappropriate smoothing length
         */
        if (fabsf(h_new - h_old) > eps * h_old) {

          /* Ok, correct then */

          /* Case where we have been oscillating around the solution */
          if ((h_new == left[i] && h_old == right[i]) ||
              (h_old == left[i] && h_new == right[i])) {

            /* Bissect the remaining interval */
            bp->h = pow_inv_dimension(
                0.5f * (pow_dimension(left[i]) + pow_dimension(right[i])));

          } else {

            /* Normal case */
            bp->h = h_new;
          }

          /* If below the absolute maximum, try again */
          if (bp->h < black_holes_h_max && bp->h > black_holes_h_min) {

            /* Flag for another round of fun */
            sid[redo] = sid[i];
            h_0[redo] = h_0[i];
            left[redo] = left[i];
            right[redo] = right[i];
            redo += 1;

            /* Re-initialise everything */
            black_holes_init_bpart(bp);

            /* Off we go ! */
            continue;

          } else if (bp->h <= black_holes_h_min) {

            /* Ok, this particle is a lost cause... */
            bp->h = black_holes_h_min;

          } else if (bp->h >= black_holes_h_max) {

            /* Ok, this particle is a lost cause... */
            bp->h = black_holes_h_max;

            /* Do some damage control if no neighbours at all were found */
            if (has_no_neighbours) {
              black_holes_bpart_has_no_neighbours(bp, cosmo);
            }

          } else {
            error(
                "Fundamental problem with the smoothing length iteration "
                "logic.");
          }
        }

        /* We now have a particle whose smoothing length has converged */

        black_holes_reset_feedback(bp);

        /* Check if h_max has increased */
        h_max = max(h_max, bp->h);
      }

      /* We now need to treat the particles whose smoothing length had not
       * converged again */

      /* Re-set the counter for the next loop (potentially). */
      bcount = redo;
      if (bcount > 0) {

        /* Climb up the cell hierarchy. */
        for (struct cell *finger = c; finger != NULL; finger = finger->parent) {

          /* Run through this cell's density interactions. */
          for (struct link *l = finger->black_holes.density; l != NULL;
               l = l->next) {

#ifdef SWIFT_DEBUG_CHECKS
            if (l->t->ti_run < r->e->ti_current)
              error("Density task should have been run.");
#endif

            /* Self-interaction? */
            if (l->t->type == task_type_self)
              runner_doself_subset_branch_bh_density(r, finger, bparts, sid,
                                                     bcount);

            /* Otherwise, pair interaction? */
            else if (l->t->type == task_type_pair) {

              /* Left or right? */
              if (l->t->ci == finger)
                runner_dopair_subset_branch_bh_density(r, finger, bparts, sid,
                                                       bcount, l->t->cj);
              else
                runner_dopair_subset_branch_bh_density(r, finger, bparts, sid,
                                                       bcount, l->t->ci);
            }

            /* Otherwise, sub-self interaction? */
            else if (l->t->type == task_type_sub_self)
              runner_dosub_subset_bh_density(r, finger, bparts, sid, bcount,
                                             NULL, 1);

            /* Otherwise, sub-pair interaction? */
            else if (l->t->type == task_type_sub_pair) {

              /* Left or right? */
              if (l->t->ci == finger)
                runner_dosub_subset_bh_density(r, finger, bparts, sid, bcount,
                                               l->t->cj, 1);
              else
                runner_dosub_subset_bh_density(r, finger, bparts, sid, bcount,
                                               l->t->ci, 1);
            }
          }
        }
      }
    }

    if (bcount) {
      error("Smoothing length failed to converge on %i particles.", bcount);
    }

    /* Be clean */
    free(left);
    free(right);
    free(sid);
    free(h_0);
  }

  /* Update h_max */
  c->black_holes.h_max = h_max;

  /* The ghost may not always be at the top level.
   * Therefore we need to update h_max between the super- and top-levels */
  if (c->black_holes.density_ghost) {
    for (struct cell *tmp = c->parent; tmp != NULL; tmp = tmp->parent) {
      atomic_max_f(&tmp->black_holes.h_max, h_max);
    }
  }

  if (timer) TIMER_TOC(timer_do_black_holes_ghost);
}

/**
 * @brief Intermediate task after the BHs have done their swallowing step.
 * This is used to update the BH quantities if necessary.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_black_holes_swallow_ghost(struct runner *r, struct cell *c,
                                         int timer) {

  struct bpart *restrict bparts = c->black_holes.parts;
  const int count = c->black_holes.count;
  const struct engine *e = r->e;
  const int with_cosmology = e->policy & engine_policy_cosmology;

  TIMER_TIC;

  /* Anything to do here? */
  if (c->black_holes.count == 0) return;
  if (!cell_is_active_black_holes(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL)
        runner_do_black_holes_swallow_ghost(r, c->progeny[k], 0);
  } else {

    /* Loop over the parts in this cell. */
    for (int i = 0; i < count; i++) {

      /* Get a direct pointer on the part. */
      struct bpart *bp = &bparts[i];

      if (bpart_is_active(bp, e)) {

        /* Get particle time-step */
        double dt;
        if (with_cosmology) {
          const integertime_t ti_step = get_integer_timestep(bp->time_bin);
          const integertime_t ti_begin =
              get_integer_time_begin(e->ti_current - 1, bp->time_bin);

          dt = cosmology_get_delta_time(e->cosmology, ti_begin,
                                        ti_begin + ti_step);
        } else {
          dt = get_timestep(bp->time_bin, e->time_base);
        }

        /* Compute the final operations for repositioning of this BH */
        black_holes_end_reposition(bp, e->black_holes_properties,
                                   e->physical_constants, e->cosmology, dt);

        /* Compute variables required for the feedback loop */
        black_holes_prepare_feedback(
            bp, e->black_holes_properties, e->physical_constants, e->cosmology,
            e->cooling_func, e->entropy_floor, e->time, with_cosmology, dt);
      }
    }
  }

  if (timer) TIMER_TOC(timer_do_black_holes_ghost);
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
  const integertime_t ti_current = e->ti_current;
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
          const integertime_t ti_begin =
              get_integer_time_begin(ti_current - 1, p->time_bin);

          dt_alpha =
              cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
        } else {
          dt_alpha = get_timestep(p->time_bin, time_base);
        }

        /* Compute variables required for the force loop */
        hydro_prepare_force(p, xp, cosmo, hydro_props, dt_alpha);
        timestep_limiter_prepare_force(p, xp);

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
  const struct hydro_props *hydro_props = e->hydro_properties;

  const int with_cosmology = (e->policy & engine_policy_cosmology);

  const float hydro_h_max = e->hydro_properties->h_max;
  const float hydro_h_min = e->hydro_properties->h_min;
  const float eps = e->hydro_properties->h_tolerance;
  const float hydro_eta_dim =
      pow_dimension(e->hydro_properties->eta_neighbours);
  const int use_mass_weighted_num_ngb =
      e->hydro_properties->use_mass_weighted_num_ngb;
  const int max_smoothing_iter = e->hydro_properties->max_smoothing_iterations;
  int redo = 0, count = 0;

  /* Running value of the maximal smoothing length */
  double h_max = c->hydro.h_max;

  TIMER_TIC;

  /* Anything to do here? */
  if (c->hydro.count == 0) return;
  if (!cell_is_active_hydro(c, e)) return;

  /* Recurse? */
  if (c->split) {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        runner_do_ghost(r, c->progeny[k], 0);

        /* Update h_max */
        h_max = max(h_max, c->progeny[k]->hydro.h_max);
      }
    }
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
          pressure_floor_end_density(p, cosmo);
          star_formation_end_density(p, xp, star_formation, cosmo);

          /* Are we using the alternative definition of the
             number of neighbours? */
          if (use_mass_weighted_num_ngb) {
#if defined(GIZMO_MFV_SPH) || defined(GIZMO_MFM_SPH) || defined(SHADOWFAX_SPH)
            error(
                "Can't use alternative neighbour definition with this scheme!");
#else
            const float inv_mass = 1.f / hydro_get_mass(p);
            p->density.wcount = p->rho * inv_mass;
            p->density.wcount_dh = p->density.rho_dh * inv_mass;
#endif
          }

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
           * smaller but has already reached the minimum). So, just tidy up
           * as if the smoothing length had converged correctly  */

#ifdef EXTRA_HYDRO_LOOP

            /* As of here, particle gradient variables will be set. */
            /* The force variables are set in the extra ghost. */

            /* Compute variables required for the gradient loop */
            hydro_prepare_gradient(p, xp, cosmo, hydro_props);

            /* The particle gradient values are now set.  Do _NOT_
               try to read any particle density variables! */

            /* Prepare the particle for the gradient loop over neighbours
             */
            hydro_reset_gradient(p);

#else
            /* Calculate the time-step for passing to hydro_prepare_force, used
             * for the evolution of alpha factors (i.e. those involved in the
             * artificial viscosity and thermal conduction terms) */
            const double time_base = e->time_base;
            const integertime_t ti_current = e->ti_current;
            double dt_alpha;

            if (with_cosmology) {
              const integertime_t ti_step = get_integer_timestep(p->time_bin);
              const integertime_t ti_begin =
                  get_integer_time_begin(ti_current - 1, p->time_bin);

              dt_alpha =
                  cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
            } else {
              dt_alpha = get_timestep(p->time_bin, time_base);
            }

            /* As of here, particle force variables will be set. */

            /* Compute variables required for the force loop */
            hydro_prepare_force(p, xp, cosmo, hydro_props, dt_alpha);
            timestep_limiter_prepare_force(p, xp);

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

        /* Check whether the particle has an inappropriate smoothing length
         */
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
            pressure_floor_init_part(p, xp);
            star_formation_init_part(p, star_formation);
            tracers_after_init(p, xp, e->internal_units, e->physical_constants,
                               with_cosmology, e->cosmology,
                               e->hydro_properties, e->cooling_func, e->time);

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
              pressure_floor_part_has_no_neighbours(p, xp, cosmo);
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

        /* Check if h_max is increased */
        h_max = max(h_max, p->h);

#ifdef EXTRA_HYDRO_LOOP

        /* As of here, particle gradient variables will be set. */
        /* The force variables are set in the extra ghost. */

        /* Compute variables required for the gradient loop */
        hydro_prepare_gradient(p, xp, cosmo, hydro_props);

        /* The particle gradient values are now set.  Do _NOT_
           try to read any particle density variables! */

        /* Prepare the particle for the gradient loop over neighbours */
        hydro_reset_gradient(p);

#else

        /* Calculate the time-step for passing to hydro_prepare_force, used
         * for the evolution of alpha factors (i.e. those involved in the
         * artificial viscosity and thermal conduction terms) */
        const double time_base = e->time_base;
        const integertime_t ti_current = e->ti_current;
        double dt_alpha;

        if (with_cosmology) {
          const integertime_t ti_step = get_integer_timestep(p->time_bin);
          const integertime_t ti_begin =
              get_integer_time_begin(ti_current - 1, p->time_bin);

          dt_alpha =
              cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
        } else {
          dt_alpha = get_timestep(p->time_bin, time_base);
        }

        /* As of here, particle force variables will be set. */

        /* Compute variables required for the force loop */
        hydro_prepare_force(p, xp, cosmo, hydro_props, dt_alpha);
        timestep_limiter_prepare_force(p, xp);

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
                                          1);

            /* Otherwise, sub-pair interaction? */
            else if (l->t->type == task_type_sub_pair) {

              /* Left or right? */
              if (l->t->ci == finger)
                runner_dosub_subset_density(r, finger, parts, pid, count,
                                            l->t->cj, 1);
              else
                runner_dosub_subset_density(r, finger, parts, pid, count,
                                            l->t->ci, 1);
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

  /* Update h_max */
  c->hydro.h_max = h_max;

  /* The ghost may not always be at the top level.
   * Therefore we need to update h_max between the super- and top-levels */
  if (c->hydro.ghost) {
    for (struct cell *tmp = c->parent; tmp != NULL; tmp = tmp->parent) {
      atomic_max_f(&tmp->hydro.h_max, h_max);
    }
  }

  if (timer) TIMER_TOC(timer_do_ghost);
}

/**
 * @brief Intermediate task after the density loop to check that the smoothing
 * lengths for dark matter particles are correct.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_dark_matter_density_ghost(struct runner *r, struct cell *c) {
    
    struct dmpart *restrict dmparts = c->dark_matter.parts;
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    const int with_cosmology = (e->policy & engine_policy_cosmology);
    
    const float dark_matter_h_max = e->sidm_properties->h_max;
    const float dark_matter_h_min = e->sidm_properties->h_min;
    const float eps = e->sidm_properties->h_tolerance;
    const float dark_matter_eta_dim =
    pow_dimension(e->sidm_properties->eta_neighbours);
    const int use_mass_weighted_num_ngb =
    e->sidm_properties->use_mass_weighted_num_ngb;
    const int max_smoothing_iter = e->sidm_properties->max_smoothing_iterations;
    int redo = 0, count = 0;
    
    /* Running value of the maximal smoothing length */
    double h_max = c->dark_matter.h_max;
    
    TIMER_TIC;
    message("Runner DM density ghost");
    
    /* Anything to do here? */
    if (c->dark_matter.count == 0) return;
    if (!cell_is_active_dark_matter(c, e)) return;
    
    /* Recurse? */
    if (c->split) {
        for (int k = 0; k < 8; k++) {
            if (c->progeny[k] != NULL) {
                runner_do_dark_matter_density_ghost(r, c->progeny[k]);
                
                /* Update h_max */
                h_max = max(h_max, c->progeny[k]->dark_matter.h_max);
            }
        }
    } else {
        
        /* Init the list of active particles that have to be updated and their
         * current smoothing lengths. */
        int *pid = NULL;
        float *h_0 = NULL;
        float *left = NULL;
        float *right = NULL;
        if ((pid = (int *)malloc(sizeof(int) * c->dark_matter.count)) == NULL)
            error("Can't allocate memory for pid.");
        if ((h_0 = (float *)malloc(sizeof(float) * c->dark_matter.count)) == NULL)
            error("Can't allocate memory for h_0.");
        if ((left = (float *)malloc(sizeof(float) * c->dark_matter.count)) == NULL)
            error("Can't allocate memory for left.");
        if ((right = (float *)malloc(sizeof(float) * c->dark_matter.count)) == NULL)
            error("Can't allocate memory for right.");
        for (int k = 0; k < c->dark_matter.count; k++)
            if (dmpart_is_active(&dmparts[k], e)) {
                pid[count] = k;
                h_0[count] = dmparts[k].h;
                left[count] = 0.f;
                right[count] = dark_matter_h_max;
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
                struct dmpart *p = &dmparts[pid[i]];
                
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
                    dark_matter_end_density(p, cosmo);
                    
                    /* Are we using the alternative definition of the
                     number of neighbours? */
                    if (use_mass_weighted_num_ngb) {
                        const float inv_mass = 1.f / p->mass;
                        p->density.wcount = p->rho * inv_mass;
                    }
                    
                    /* Compute one step of the Newton-Raphson scheme */
                    const float n_sum = p->density.wcount * h_old_dim;
                    const float n_target = dark_matter_eta_dim;
                    const float f = n_sum - n_target;
                    const float f_prime = p->density.wcount_dh * h_old_dim +
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
                    if (((p->h >= dark_matter_h_max) && (f < 0.f)) ||
                        ((p->h <= dark_matter_h_min) && (f > 0.f))) {
                        
                        /* We have a particle whose smoothing length is already set (wants
                         * to be larger but has already hit the maximum OR wants to be
                         * smaller but has already reached the minimum). So, just tidy up
                         * as if the smoothing length had converged correctly  */
                        

                        /* Calculate the time-step for passing to hydro_prepare_force, used
                         * for the evolution of alpha factors (i.e. those involved in the
                         * artificial viscosity and thermal conduction terms) */
                        const double time_base = e->time_base;
                        const integertime_t ti_current = e->ti_current;
                        double dt_alpha;
                        
                        if (with_cosmology) {
                            const integertime_t ti_step = get_integer_timestep(p->time_bin);
                            const integertime_t ti_begin =
                            get_integer_time_begin(ti_current - 1, p->time_bin);
                            
                            dt_alpha =
                            cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
                        } else {
                            dt_alpha = get_timestep(p->time_bin, time_base);
                        }
                        
                        /* As of here, particle force variables will be set. */
                        
                        /* Compute variables required for the force loop */
                        /*hydro_prepare_force(p, xp, cosmo, hydro_props, dt_alpha);
                        timestep_limiter_prepare_force(p, xp);*/
                        
                        /* The particle force values are now set.  Do _NOT_
                         try to read any particle density variables! */
                        
                        /* Prepare the particle for the force loop over neighbours */
                        /*hydro_reset_acceleration(p);*/
                        
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
                                num_reruns, p->id_or_neg_offset, h_init, h_old, h_new, f, f_prime, n_sum,
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
                
                /* Check whether the particle has an inappropriate smoothing length
                 */
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
                    if (p->h < dark_matter_h_max && p->h > dark_matter_h_min) {
                        
                        /* Flag for another round of fun */
                        pid[redo] = pid[i];
                        h_0[redo] = h_0[i];
                        left[redo] = left[i];
                        right[redo] = right[i];
                        redo += 1;
                        
                        /* Re-initialise everything */
                        dark_matter_init_dmpart(p);
                        
                        /* Off we go ! */
                        continue;
                        
                    } else if (p->h <= dark_matter_h_min) {
                        
                        /* Ok, this particle is a lost cause... */
                        p->h = dark_matter_h_min;
                        
                    } else if (p->h >= dark_matter_h_max) {
                        
                        /* Ok, this particle is a lost cause... */
                        p->h = dark_matter_h_max;
                        
                        /* Do some damage control if no neighbours at all were found */
                        if (has_no_neighbours) {
                            dark_matter_part_has_no_neighbours(p, cosmo);
                        }
                        
                    } else {
                        error(
                              "Fundamental problem with the smoothing length iteration "
                              "logic.");
                    }
                }
                
                /* We now have a particle whose smoothing length has converged */
                
                /* Check if h_max is increased */
                h_max = max(h_max, p->h);
                
                
                /* Calculate the time-step for passing to hydro_prepare_force, used
                 * for the evolution of alpha factors (i.e. those involved in the
                 * artificial viscosity and thermal conduction terms) */
                const double time_base = e->time_base;
                const integertime_t ti_current = e->ti_current;
                double dt_alpha;
                
                if (with_cosmology) {
                    const integertime_t ti_step = get_integer_timestep(p->time_bin);
                    const integertime_t ti_begin =
                    get_integer_time_begin(ti_current - 1, p->time_bin);
                    
                    dt_alpha =
                    cosmology_get_delta_time(cosmo, ti_begin, ti_begin + ti_step);
                } else {
                    dt_alpha = get_timestep(p->time_bin, time_base);
                }
                
                /* As of here, particle force variables will be set. */
                
                /* Compute variables required for the force loop */
                /*hydro_prepare_force(p, xp, cosmo, hydro_props, dt_alpha);
                timestep_limiter_prepare_force(p, xp);*/
                
                /* The particle force values are now set.  Do _NOT_
                 try to read any particle density variables! */
                
                /* Prepare the particle for the force loop over neighbours */
                /*hydro_reset_acceleration(p);*/
                
            }
            
            /* We now need to treat the particles whose smoothing length had not
             * converged again */
            
            /* Re-set the counter for the next loop (potentially). */
            count = redo;
            if (count > 0) {
                
                /* Climb up the cell hierarchy. */
                for (struct cell *finger = c; finger != NULL; finger = finger->parent) {
                    
                    /* Run through this cell's density interactions. */
                    for (struct link *l = finger->dark_matter.density; l != NULL; l = l->next) {
                        
#ifdef SWIFT_DEBUG_CHECKS
                        if (l->t->ti_run < r->e->ti_current)
                            error("Density task should have been run.");
#endif
                        
                        /* Self-interaction? */
                        if (l->t->type == task_type_self)
                            runner_doself_subset_dark_matter_density(r, finger, dmparts, pid, count);
                        
                        /* Otherwise, pair interaction? */
                        else if (l->t->type == task_type_pair) {
                            
                            /* Left or right? */
                            if (l->t->ci == finger)
                                runner_dopair_subset_dark_matter_density(r, finger, dmparts, pid, l->t->cj, count);
                            else
                                runner_dopair_subset_dark_matter_density(r, finger, dmparts, pid, l->t->ci, count);
                        }
                        
                        /* Otherwise, sub-self interaction? */
                        else if (l->t->type == task_type_sub_self)
                            runner_doself_subset_dark_matter_density(r, finger, dmparts, pid, count);
                        
                        /* Otherwise, sub-pair interaction? */
                        else if (l->t->type == task_type_sub_pair) {
                            
                            /* Left or right? */
                            if (l->t->ci == finger)
                                runner_dopair_subset_dark_matter_density(r, finger, dmparts, pid, l->t->cj, count);
                            else
                                runner_dopair_subset_dark_matter_density(r, finger, dmparts, pid, l->t->ci, count);
                        }
                    }
                }
            }
        }
        
        if (count) {
            error("Smoothing length failed to converge on %i DM particles.", count);
        }
        
        /* Be clean */
        free(left);
        free(right);
        free(pid);
        free(h_0);
    }
    
    /* Update h_max */
    c->dark_matter.h_max = h_max;
    
    /* The ghost may not always be at the top level.
     * Therefore we need to update h_max between the super- and top-levels */
    if (c->dark_matter.ghost) {
        for (struct cell *tmp = c->parent; tmp != NULL; tmp = tmp->parent) {
            atomic_max_f(&tmp->dark_matter.h_max, h_max);
        }
    }
}
