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
#include "cell.h"

/* Local headers. */
#include "active.h"
#include "drift.h"
#include "feedback.h"
#include "gravity.h"
#include "multipole.h"
#include "pressure_floor.h"
#include "rt.h"
#include "star_formation.h"
#include "tracers.h"

/**
 * @brief Recursively drifts the #part in a cell hierarchy.
 *
 * @param c The #cell.
 * @param e The #engine (to get ti_current).
 * @param force Drift the particles irrespective of the #cell flags.
 */
void cell_drift_part(struct cell *c, const struct engine *e, int force) {
  const int periodic = e->s->periodic;
  const double dim[3] = {e->s->dim[0], e->s->dim[1], e->s->dim[2]};
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const float hydro_h_max = e->hydro_properties->h_max;
  const float hydro_h_min = e->hydro_properties->h_min;
  const integertime_t ti_old_part = c->hydro.ti_old_part;
  const integertime_t ti_current = e->ti_current;
  struct part *const parts = c->hydro.parts;
  struct xpart *const xparts = c->hydro.xparts;

  float dx_max = 0.f, dx2_max = 0.f;
  float dx_max_sort = 0.0f, dx2_max_sort = 0.f;
  float cell_h_max = 0.f;

  /* Drift irrespective of cell flags? */
  force = (force || cell_get_flag(c, cell_flag_do_hydro_drift));

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we only drift local cells. */
  if (c->nodeID != engine_rank) error("Drifting a foreign cell is nope.");

  /* Check that we are actually going to move forward. */
  if (ti_current < ti_old_part) error("Attempt to drift to the past");
#endif

  /* Early abort? */
  if (c->hydro.count == 0) {
    /* Clear the drift flags. */
    cell_clear_flag(c, cell_flag_do_hydro_drift | cell_flag_do_hydro_sub_drift);

    /* Update the time of the last drift */
    c->hydro.ti_old_part = ti_current;

    return;
  }

  /* Ok, we have some particles somewhere in the hierarchy to drift */

  /* Are we not in a leaf ? */
  if (c->split && (force || cell_get_flag(c, cell_flag_do_hydro_sub_drift))) {

    /* Loop over the progeny and collect their data. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];

        /* Collect */
        cell_drift_part(cp, e, force);

        /* Update */
        dx_max = max(dx_max, cp->hydro.dx_max_part);
        dx_max_sort = max(dx_max_sort, cp->hydro.dx_max_sort);
        cell_h_max = max(cell_h_max, cp->hydro.h_max);
      }
    }

    /* Store the values */
    c->hydro.h_max = cell_h_max;
    c->hydro.dx_max_part = dx_max;
    c->hydro.dx_max_sort = dx_max_sort;

    /* Update the time of the last drift */
    c->hydro.ti_old_part = ti_current;

  } else if (!c->split && force && ti_current > ti_old_part) {
    /* Drift from the last time the cell was drifted to the current time */
    double dt_drift, dt_kick_grav, dt_kick_hydro, dt_therm;
    if (with_cosmology) {
      dt_drift =
          cosmology_get_drift_factor(e->cosmology, ti_old_part, ti_current);
      dt_kick_grav =
          cosmology_get_grav_kick_factor(e->cosmology, ti_old_part, ti_current);
      dt_kick_hydro = cosmology_get_hydro_kick_factor(e->cosmology, ti_old_part,
                                                      ti_current);
      dt_therm = cosmology_get_therm_kick_factor(e->cosmology, ti_old_part,
                                                 ti_current);
    } else {
      dt_drift = (ti_current - ti_old_part) * e->time_base;
      dt_kick_grav = (ti_current - ti_old_part) * e->time_base;
      dt_kick_hydro = (ti_current - ti_old_part) * e->time_base;
      dt_therm = (ti_current - ti_old_part) * e->time_base;
    }

    /* Loop over all the gas particles in the cell */
    const size_t nr_parts = c->hydro.count;
    for (size_t k = 0; k < nr_parts; k++) {
      /* Get a handle on the part. */
      struct part *const p = &parts[k];
      struct xpart *const xp = &xparts[k];

      /* Ignore inhibited particles */
      if (part_is_inhibited(p, e)) continue;

      /* Apply the effects of feedback on this particle
       * (Note: Only used in schemes that have a delayed feedback mechanism
       * otherwise just an empty function) */
      feedback_update_part(p, xp, e);

      /* Drift... */
      drift_part(p, xp, dt_drift, dt_kick_hydro, dt_kick_grav, dt_therm,
                 ti_old_part, ti_current, e->cosmology, e->hydro_properties,
                 e->entropy_floor);

      /* Update the tracers properties */
      tracers_after_drift(p, xp, e->internal_units, e->physical_constants,
                          with_cosmology, e->cosmology, e->hydro_properties,
                          e->cooling_func, e->time);

#ifdef SWIFT_DEBUG_CHECKS
      /* Make sure the particle does not drift by more than a box length. */
      if (fabs(xp->v_full[0] * dt_drift) > e->s->dim[0] ||
          fabs(xp->v_full[1] * dt_drift) > e->s->dim[1] ||
          fabs(xp->v_full[2] * dt_drift) > e->s->dim[2]) {
        error(
            "Particle drifts by more than a box length! id %llu xp->v_full "
            "%.5e %.5e %.5e p->v %.5e %.5e %.5e",
            p->id, xp->v_full[0], xp->v_full[1], xp->v_full[2], p->v[0],
            p->v[1], p->v[2]);
      }
#endif

      /* In non-periodic BC runs, remove particles that crossed the border */
      if (!periodic) {

        /* Did the particle leave the box?  */
        if ((p->x[0] > dim[0]) || (p->x[0] < 0.) ||  // x
            (p->x[1] > dim[1]) || (p->x[1] < 0.) ||  // y
            (p->x[2] > dim[2]) || (p->x[2] < 0.)) {  // z

          lock_lock(&e->s->lock);

          /* Re-check that the particle has not been removed
           * by another thread before we do the deed. */
          if (!part_is_inhibited(p, e)) {

#ifdef WITH_LOGGER
            if (e->policy & engine_policy_logger) {
              /* Log the particle one last time. */
              logger_log_part(
                  e->logger, p, xp, e, /* log_all */ 1,
                  logger_pack_flags_and_data(logger_flag_delete, 0));
            }
#endif

            /* One last action before death? */
            hydro_remove_part(p, xp);

            /* Remove the particle entirely */
            cell_remove_part(e, c, p, xp);
          }

          if (lock_unlock(&e->s->lock) != 0)
            error("Failed to unlock the space!");

          continue;
        }
      }

      /* Limit h to within the allowed range */
      p->h = min(p->h, hydro_h_max);
      p->h = max(p->h, hydro_h_min);

      /* Compute (square of) motion since last cell construction */
      const float dx2 = xp->x_diff[0] * xp->x_diff[0] +
                        xp->x_diff[1] * xp->x_diff[1] +
                        xp->x_diff[2] * xp->x_diff[2];
      dx2_max = max(dx2_max, dx2);
      const float dx2_sort = xp->x_diff_sort[0] * xp->x_diff_sort[0] +
                             xp->x_diff_sort[1] * xp->x_diff_sort[1] +
                             xp->x_diff_sort[2] * xp->x_diff_sort[2];
      dx2_max_sort = max(dx2_max_sort, dx2_sort);

      /* Update the maximal smoothing length in the cell */
      cell_h_max = max(cell_h_max, p->h);

      /* Mark the particle has not being swallowed */
      black_holes_mark_part_as_not_swallowed(&p->black_holes_data);

      /* Get ready for a density calculation */
      if (part_is_active(p, e)) {
        hydro_init_part(p, &e->s->hs);
        black_holes_init_potential(&p->black_holes_data);
        chemistry_init_part(p, e->chemistry);
        pressure_floor_init_part(p, xp);
        star_formation_init_part(p, e->star_formation);
        tracers_after_init(p, xp, e->internal_units, e->physical_constants,
                           with_cosmology, e->cosmology, e->hydro_properties,
                           e->cooling_func, e->time);
        rt_init_part(p);
        rt_reset_part(p); /* also reset density unrelated RT part quantities */
      }
    }

    /* Now, get the maximal particle motion from its square */
    dx_max = sqrtf(dx2_max);
    dx_max_sort = sqrtf(dx2_max_sort);

    /* Store the values */
    c->hydro.h_max = cell_h_max;
    c->hydro.dx_max_part = dx_max;
    c->hydro.dx_max_sort = dx_max_sort;

    /* Update the time of the last drift */
    c->hydro.ti_old_part = ti_current;
  }

  /* Clear the drift flags. */
  cell_clear_flag(c, cell_flag_do_hydro_drift | cell_flag_do_hydro_sub_drift);
}

/**
 * @brief Recursively drifts the #gpart in a cell hierarchy.
 *
 * @param c The #cell.
 * @param e The #engine (to get ti_current).
 * @param force Drift the particles irrespective of the #cell flags.
 */
void cell_drift_gpart(struct cell *c, const struct engine *e, int force) {
  const int periodic = e->s->periodic;
  const double dim[3] = {e->s->dim[0], e->s->dim[1], e->s->dim[2]};
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const integertime_t ti_old_gpart = c->grav.ti_old_part;
  const integertime_t ti_current = e->ti_current;
  struct gpart *const gparts = c->grav.parts;
  const struct gravity_props *grav_props = e->gravity_properties;

  /* Drift irrespective of cell flags? */
  force = (force || cell_get_flag(c, cell_flag_do_grav_drift));

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we only drift local cells. */
  if (c->nodeID != engine_rank) error("Drifting a foreign cell is nope.");

  /* Check that we are actually going to move forward. */
  if (ti_current < ti_old_gpart) error("Attempt to drift to the past");
#endif

  /* Early abort? */
  if (c->grav.count == 0) {
    /* Clear the drift flags. */
    cell_clear_flag(c, cell_flag_do_grav_drift | cell_flag_do_grav_sub_drift);

    /* Update the time of the last drift */
    c->grav.ti_old_part = ti_current;

    return;
  }

  /* Ok, we have some particles somewhere in the hierarchy to drift */

  /* Are we not in a leaf ? */
  if (c->split && (force || cell_get_flag(c, cell_flag_do_grav_sub_drift))) {

    /* Loop over the progeny and collect their data. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];

        /* Recurse */
        cell_drift_gpart(cp, e, force);
      }
    }

    /* Update the time of the last drift */
    c->grav.ti_old_part = ti_current;

  } else if (!c->split && force && ti_current > ti_old_gpart) {
    /* Drift from the last time the cell was drifted to the current time */
    double dt_drift;
    if (with_cosmology) {
      dt_drift =
          cosmology_get_drift_factor(e->cosmology, ti_old_gpart, ti_current);
    } else {
      dt_drift = (ti_current - ti_old_gpart) * e->time_base;
    }

    /* Loop over all the g-particles in the cell */
    const size_t nr_gparts = c->grav.count;
    for (size_t k = 0; k < nr_gparts; k++) {
      /* Get a handle on the gpart. */
      struct gpart *const gp = &gparts[k];

      /* Ignore inhibited particles */
      if (gpart_is_inhibited(gp, e)) continue;

      /* Drift... */
      drift_gpart(gp, dt_drift, ti_old_gpart, ti_current, grav_props, e);

#ifdef SWIFT_DEBUG_CHECKS
      /* Make sure the particle does not drift by more than a box length. */
      if (fabs(gp->v_full[0] * dt_drift) > e->s->dim[0] ||
          fabs(gp->v_full[1] * dt_drift) > e->s->dim[1] ||
          fabs(gp->v_full[2] * dt_drift) > e->s->dim[2]) {
        error(
            "Particle drifts by more than a box length! gp->v_full %.5e %.5e "
            "%.5e",
            gp->v_full[0], gp->v_full[1], gp->v_full[2]);
      }
#endif

      /* In non-periodic BC runs, remove particles that crossed the border */
      if (!periodic) {

        /* Did the particle leave the box?  */
        if ((gp->x[0] > dim[0]) || (gp->x[0] < 0.) ||  // x
            (gp->x[1] > dim[1]) || (gp->x[1] < 0.) ||  // y
            (gp->x[2] > dim[2]) || (gp->x[2] < 0.)) {  // z

          lock_lock(&e->s->lock);

          /* Re-check that the particle has not been removed
           * by another thread before we do the deed. */
          if (!gpart_is_inhibited(gp, e)) {

            /* Remove the particle entirely */
            if (gp->type == swift_type_dark_matter) {

#ifdef WITH_LOGGER
              if (e->policy & engine_policy_logger) {
                /* Log the particle one last time. */
                logger_log_gpart(
                    e->logger, gp, e, /* log_all */ 1,
                    logger_pack_flags_and_data(logger_flag_delete, 0));
              }
#endif

              /* Remove the particle */
              cell_remove_gpart(e, c, gp);
            }
          }

          if (lock_unlock(&e->s->lock) != 0)
            error("Failed to unlock the space!");

          continue;
        }
      }

      /* Init gravity force fields. */
      if (gpart_is_active(gp, e)) {
        gravity_init_gpart(gp);
      }
    }

    /* Update the time of the last drift */
    c->grav.ti_old_part = ti_current;
  }

  /* Clear the drift flags. */
  cell_clear_flag(c, cell_flag_do_grav_drift | cell_flag_do_grav_sub_drift);
}

/**
 * @brief Recursively drifts the #spart in a cell hierarchy.
 *
 * @param c The #cell.
 * @param e The #engine (to get ti_current).
 * @param force Drift the particles irrespective of the #cell flags.
 */
void cell_drift_spart(struct cell *c, const struct engine *e, int force) {
  const int periodic = e->s->periodic;
  const double dim[3] = {e->s->dim[0], e->s->dim[1], e->s->dim[2]};
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const int with_rt = (e->policy & engine_policy_rt);
  const float stars_h_max = e->hydro_properties->h_max;
  const float stars_h_min = e->hydro_properties->h_min;
  const integertime_t ti_old_spart = c->stars.ti_old_part;
  const integertime_t ti_current = e->ti_current;
  struct spart *const sparts = c->stars.parts;

  float dx_max = 0.f, dx2_max = 0.f;
  float dx_max_sort = 0.0f, dx2_max_sort = 0.f;
  float cell_h_max = 0.f;

  /* Drift irrespective of cell flags? */
  force = (force || cell_get_flag(c, cell_flag_do_stars_drift));

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we only drift local cells. */
  if (c->nodeID != engine_rank) error("Drifting a foreign cell is nope.");

  /* Check that we are actually going to move forward. */
  if (ti_current < ti_old_spart) error("Attempt to drift to the past");
#endif

  /* Early abort? */
  if (c->stars.count == 0) {
    /* Clear the drift flags. */
    cell_clear_flag(c, cell_flag_do_stars_drift | cell_flag_do_stars_sub_drift);

    /* Update the time of the last drift */
    c->stars.ti_old_part = ti_current;

    return;
  }

  /* Ok, we have some particles somewhere in the hierarchy to drift */

  /* Are we not in a leaf ? */
  if (c->split && (force || cell_get_flag(c, cell_flag_do_stars_sub_drift))) {

    /* Loop over the progeny and collect their data. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];

        /* Recurse */
        cell_drift_spart(cp, e, force);

        /* Update */
        dx_max = max(dx_max, cp->stars.dx_max_part);
        dx_max_sort = max(dx_max_sort, cp->stars.dx_max_sort);
        cell_h_max = max(cell_h_max, cp->stars.h_max);
      }
    }

    /* Store the values */
    c->stars.h_max = cell_h_max;
    c->stars.dx_max_part = dx_max;
    c->stars.dx_max_sort = dx_max_sort;

    /* Update the time of the last drift */
    c->stars.ti_old_part = ti_current;

  } else if (!c->split && force && ti_current > ti_old_spart) {
    /* Drift from the last time the cell was drifted to the current time */
    double dt_drift;
    if (with_cosmology) {
      dt_drift =
          cosmology_get_drift_factor(e->cosmology, ti_old_spart, ti_current);
    } else {
      dt_drift = (ti_current - ti_old_spart) * e->time_base;
    }

    /* Loop over all the star particles in the cell */
    const size_t nr_sparts = c->stars.count;
    for (size_t k = 0; k < nr_sparts; k++) {
      /* Get a handle on the spart. */
      struct spart *const sp = &sparts[k];

      /* Ignore inhibited particles */
      if (spart_is_inhibited(sp, e)) continue;

      /* Drift... */
      drift_spart(sp, dt_drift, ti_old_spart, ti_current);

#ifdef SWIFT_DEBUG_CHECKS
      /* Make sure the particle does not drift by more than a box length. */
      if (fabs(sp->v[0] * dt_drift) > e->s->dim[0] ||
          fabs(sp->v[1] * dt_drift) > e->s->dim[1] ||
          fabs(sp->v[2] * dt_drift) > e->s->dim[2]) {
        error("Particle drifts by more than a box length!");
      }
#endif

      /* In non-periodic BC runs, remove particles that crossed the border */
      if (!periodic) {

        /* Did the particle leave the box?  */
        if ((sp->x[0] > dim[0]) || (sp->x[0] < 0.) ||  // x
            (sp->x[1] > dim[1]) || (sp->x[1] < 0.) ||  // y
            (sp->x[2] > dim[2]) || (sp->x[2] < 0.)) {  // z

          lock_lock(&e->s->lock);

          /* Re-check that the particle has not been removed
           * by another thread before we do the deed. */
          if (!spart_is_inhibited(sp, e)) {

#ifdef WITH_LOGGER
            if (e->policy & engine_policy_logger) {
              /* Log the particle one last time. */
              logger_log_spart(
                  e->logger, sp, e, /* log_all */ 1,
                  logger_pack_flags_and_data(logger_flag_delete, 0));
            }
#endif

            /* Remove the particle entirely */
            cell_remove_spart(e, c, sp);
          }

          if (lock_unlock(&e->s->lock) != 0)
            error("Failed to unlock the space!");

          continue;
        }
      }

      /* Limit h to within the allowed range */
      sp->h = min(sp->h, stars_h_max);
      sp->h = max(sp->h, stars_h_min);

      /* Compute (square of) motion since last cell construction */
      const float dx2 = sp->x_diff[0] * sp->x_diff[0] +
                        sp->x_diff[1] * sp->x_diff[1] +
                        sp->x_diff[2] * sp->x_diff[2];
      dx2_max = max(dx2_max, dx2);

      const float dx2_sort = sp->x_diff_sort[0] * sp->x_diff_sort[0] +
                             sp->x_diff_sort[1] * sp->x_diff_sort[1] +
                             sp->x_diff_sort[2] * sp->x_diff_sort[2];

      dx2_max_sort = max(dx2_max_sort, dx2_sort);

      /* Maximal smoothing length */
      cell_h_max = max(cell_h_max, sp->h);

      /* Get ready for a density calculation */
      if (spart_is_active(sp, e)) {
        stars_init_spart(sp);
        feedback_init_spart(sp);

        if (with_rt) {
          rt_init_spart(sp);

          /* get star's age and time step for stellar emission rates */
          const integertime_t ti_begin =
              get_integer_time_begin(ti_current - 1, sp->time_bin);
          const integertime_t ti_step = get_integer_timestep(sp->time_bin);

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
                e->cosmology, (double)sp->birth_scale_factor, e->cosmology->a);
          } else {
            star_age_end_of_step = e->time - (double)sp->birth_time;
          }

          rt_compute_stellar_emission_rate(sp, e->time, star_age_end_of_step,
                                           dt_star);
        }
      }
    }

    /* Now, get the maximal particle motion from its square */
    dx_max = sqrtf(dx2_max);
    dx_max_sort = sqrtf(dx2_max_sort);

    /* Store the values */
    c->stars.h_max = cell_h_max;
    c->stars.dx_max_part = dx_max;
    c->stars.dx_max_sort = dx_max_sort;

    /* Update the time of the last drift */
    c->stars.ti_old_part = ti_current;
  }

  /* Clear the drift flags. */
  cell_clear_flag(c, cell_flag_do_stars_drift | cell_flag_do_stars_sub_drift);
}

/**
 * @brief Recursively drifts the #bpart in a cell hierarchy.
 *
 * @param c The #cell.
 * @param e The #engine (to get ti_current).
 * @param force Drift the particles irrespective of the #cell flags.
 */
void cell_drift_bpart(struct cell *c, const struct engine *e, int force) {

  const int periodic = e->s->periodic;
  const double dim[3] = {e->s->dim[0], e->s->dim[1], e->s->dim[2]};
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const float black_holes_h_max = e->hydro_properties->h_max;
  const float black_holes_h_min = e->hydro_properties->h_min;
  const integertime_t ti_old_bpart = c->black_holes.ti_old_part;
  const integertime_t ti_current = e->ti_current;
  struct bpart *const bparts = c->black_holes.parts;

  float dx_max = 0.f, dx2_max = 0.f;
  float cell_h_max = 0.f;

  /* Drift irrespective of cell flags? */
  force = (force || cell_get_flag(c, cell_flag_do_bh_drift));

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we only drift local cells. */
  if (c->nodeID != engine_rank) error("Drifting a foreign cell is nope.");

  /* Check that we are actually going to move forward. */
  if (ti_current < ti_old_bpart) error("Attempt to drift to the past");
#endif

  /* Early abort? */
  if (c->black_holes.count == 0) {

    /* Clear the drift flags. */
    cell_clear_flag(c, cell_flag_do_bh_drift | cell_flag_do_bh_sub_drift);

    /* Update the time of the last drift */
    c->black_holes.ti_old_part = ti_current;

    return;
  }

  /* Ok, we have some particles somewhere in the hierarchy to drift */

  /* Are we not in a leaf ? */
  if (c->split && (force || cell_get_flag(c, cell_flag_do_bh_sub_drift))) {

    /* Loop over the progeny and collect their data. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];

        /* Recurse */
        cell_drift_bpart(cp, e, force);

        /* Update */
        dx_max = max(dx_max, cp->black_holes.dx_max_part);
        cell_h_max = max(cell_h_max, cp->black_holes.h_max);
      }
    }

    /* Store the values */
    c->black_holes.h_max = cell_h_max;
    c->black_holes.dx_max_part = dx_max;

    /* Update the time of the last drift */
    c->black_holes.ti_old_part = ti_current;

  } else if (!c->split && force && ti_current > ti_old_bpart) {

    /* Drift from the last time the cell was drifted to the current time */
    double dt_drift;
    if (with_cosmology) {
      dt_drift =
          cosmology_get_drift_factor(e->cosmology, ti_old_bpart, ti_current);
    } else {
      dt_drift = (ti_current - ti_old_bpart) * e->time_base;
    }

    /* Loop over all the star particles in the cell */
    const size_t nr_bparts = c->black_holes.count;
    for (size_t k = 0; k < nr_bparts; k++) {

      /* Get a handle on the bpart. */
      struct bpart *const bp = &bparts[k];

      /* Ignore inhibited particles */
      if (bpart_is_inhibited(bp, e)) continue;

      /* Drift... */
      drift_bpart(bp, dt_drift, ti_old_bpart, ti_current);

#ifdef SWIFT_DEBUG_CHECKS
      /* Make sure the particle does not drift by more than a box length. */
      if (fabs(bp->v[0] * dt_drift) > e->s->dim[0] ||
          fabs(bp->v[1] * dt_drift) > e->s->dim[1] ||
          fabs(bp->v[2] * dt_drift) > e->s->dim[2]) {
        error("Particle drifts by more than a box length!");
      }
#endif

      /* In non-periodic BC runs, remove particles that crossed the border */
      if (!periodic) {

        /* Did the particle leave the box?  */
        if ((bp->x[0] > dim[0]) || (bp->x[0] < 0.) ||  // x
            (bp->x[1] > dim[1]) || (bp->x[1] < 0.) ||  // y
            (bp->x[2] > dim[2]) || (bp->x[2] < 0.)) {  // z

          lock_lock(&e->s->lock);

          /* Re-check that the particle has not been removed
           * by another thread before we do the deed. */
          if (!bpart_is_inhibited(bp, e)) {

#ifdef WITH_LOGGER
            if (e->policy & engine_policy_logger) {
              error("Logging of black hole particles is not yet implemented.");
            }
#endif

            /* Remove the particle entirely */
            cell_remove_bpart(e, c, bp);
          }

          if (lock_unlock(&e->s->lock) != 0)
            error("Failed to unlock the space!");

          continue;
        }
      }

      /* Limit h to within the allowed range */
      bp->h = min(bp->h, black_holes_h_max);
      bp->h = max(bp->h, black_holes_h_min);

      /* Compute (square of) motion since last cell construction */
      const float dx2 = bp->x_diff[0] * bp->x_diff[0] +
                        bp->x_diff[1] * bp->x_diff[1] +
                        bp->x_diff[2] * bp->x_diff[2];
      dx2_max = max(dx2_max, dx2);

      /* Maximal smoothing length */
      cell_h_max = max(cell_h_max, bp->h);

      /* Mark the particle has not being swallowed */
      black_holes_mark_bpart_as_not_swallowed(&bp->merger_data);

      /* Get ready for a density calculation */
      if (bpart_is_active(bp, e)) {
        black_holes_init_bpart(bp);
      }
    }

    /* Now, get the maximal particle motion from its square */
    dx_max = sqrtf(dx2_max);

    /* Store the values */
    c->black_holes.h_max = cell_h_max;
    c->black_holes.dx_max_part = dx_max;

    /* Update the time of the last drift */
    c->black_holes.ti_old_part = ti_current;
  }

  /* Clear the drift flags. */
  cell_clear_flag(c, cell_flag_do_bh_drift | cell_flag_do_bh_sub_drift);
}

/**
 * @brief Recursively drifts the #sink's in a cell hierarchy.
 *
 * @param c The #cell.
 * @param e The #engine (to get ti_current).
 * @param force Drift the particles irrespective of the #cell flags.
 */
void cell_drift_sink(struct cell *c, const struct engine *e, int force) {

  const int periodic = e->s->periodic;
  const double dim[3] = {e->s->dim[0], e->s->dim[1], e->s->dim[2]};
  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const integertime_t ti_old_sink = c->sinks.ti_old_part;
  const integertime_t ti_current = e->ti_current;
  struct sink *const sinks = c->sinks.parts;

  float dx_max = 0.f, dx2_max = 0.f;
  float cell_r_max = 0.f;

  /* Drift irrespective of cell flags? */
  force = (force || cell_get_flag(c, cell_flag_do_sink_drift));

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we only drift local cells. */
  if (c->nodeID != engine_rank) error("Drifting a foreign cell is nope.");

  /* Check that we are actually going to move forward. */
  if (ti_current < ti_old_sink) error("Attempt to drift to the past");
#endif

  /* Early abort? */
  if (c->sinks.count == 0) {

    /* Clear the drift flags. */
    cell_clear_flag(c, cell_flag_do_sink_drift | cell_flag_do_sink_sub_drift);

    /* Update the time of the last drift */
    c->sinks.ti_old_part = ti_current;

    return;
  }

  /* Ok, we have some particles somewhere in the hierarchy to drift */

  /* Are we not in a leaf ? */
  if (c->split && (force || cell_get_flag(c, cell_flag_do_sink_sub_drift))) {

    /* Loop over the progeny and collect their data. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
        struct cell *cp = c->progeny[k];

        /* Recurse */
        cell_drift_sink(cp, e, force);

        /* Update */
        dx_max = max(dx_max, cp->sinks.dx_max_part);
        cell_r_max = max(cell_r_max, cp->sinks.r_cut_max);
      }
    }

    /* Store the values */
    c->sinks.r_cut_max = cell_r_max;
    c->sinks.dx_max_part = dx_max;

    /* Update the time of the last drift */
    c->sinks.ti_old_part = ti_current;

  } else if (!c->split && force && ti_current > ti_old_sink) {

    /* Drift from the last time the cell was drifted to the current time */
    double dt_drift;
    if (with_cosmology) {
      dt_drift =
          cosmology_get_drift_factor(e->cosmology, ti_old_sink, ti_current);
    } else {
      dt_drift = (ti_current - ti_old_sink) * e->time_base;
    }

    /* Loop over all the star particles in the cell */
    const size_t nr_sinks = c->sinks.count;
    for (size_t k = 0; k < nr_sinks; k++) {

      /* Get a handle on the sink. */
      struct sink *const sink = &sinks[k];

      /* Ignore inhibited particles */
      if (sink_is_inhibited(sink, e)) continue;

      /* Drift... */
      drift_sink(sink, dt_drift, ti_old_sink, ti_current);

#ifdef SWIFT_DEBUG_CHECKS
      /* Make sure the particle does not drift by more than a box length. */
      if (fabs(sink->v[0] * dt_drift) > e->s->dim[0] ||
          fabs(sink->v[1] * dt_drift) > e->s->dim[1] ||
          fabs(sink->v[2] * dt_drift) > e->s->dim[2]) {
        error("Particle drifts by more than a box length!");
      }
#endif

      /* In non-periodic BC runs, remove particles that crossed the border */
      if (!periodic) {

        /* Did the particle leave the box?  */
        if ((sink->x[0] > dim[0]) || (sink->x[0] < 0.) ||  // x
            (sink->x[1] > dim[1]) || (sink->x[1] < 0.) ||  // y
            (sink->x[2] > dim[2]) || (sink->x[2] < 0.)) {  // z

          lock_lock(&e->s->lock);

          /* Re-check that the particle has not been removed
           * by another thread before we do the deed. */
          if (!sink_is_inhibited(sink, e)) {

#ifdef WITH_LOGGER
            if (e->policy & engine_policy_logger) {
              error("Logging of sink particles is not yet implemented.");
            }
#endif

            /* Remove the particle entirely */
            // cell_remove_sink(e, c, bp);
            error("TODO: loic implement cell_remove_sink");
          }

          if (lock_unlock(&e->s->lock) != 0)
            error("Failed to unlock the space!");

          continue;
        }
      }

      /* sp->h does not need to be limited. */

      /* Compute (square of) motion since last cell construction */
      const float dx2 = sink->x_diff[0] * sink->x_diff[0] +
                        sink->x_diff[1] * sink->x_diff[1] +
                        sink->x_diff[2] * sink->x_diff[2];
      dx2_max = max(dx2_max, dx2);

      /* Maximal smoothing length */
      cell_r_max = max(cell_r_max, sink->r_cut);

      /* Get ready for a density calculation */
      if (sink_is_active(sink, e)) {
        sink_init_sink(sink);
      }
    }

    /* Now, get the maximal particle motion from its square */
    dx_max = sqrtf(dx2_max);

    /* Store the values */
    c->sinks.r_cut_max = cell_r_max;
    c->sinks.dx_max_part = dx_max;

    /* Update the time of the last drift */
    c->sinks.ti_old_part = ti_current;
  }

  /* Clear the drift flags. */
  cell_clear_flag(c, cell_flag_do_sink_drift | cell_flag_do_sink_sub_drift);
}

/**
 * @brief Recursively drifts all multipoles in a cell hierarchy.
 *
 * @param c The #cell.
 * @param e The #engine (to get ti_current).
 */
void cell_drift_all_multipoles(struct cell *c, const struct engine *e) {
  const integertime_t ti_old_multipole = c->grav.ti_old_multipole;
  const integertime_t ti_current = e->ti_current;

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we are actually going to move forward. */
  if (ti_current < ti_old_multipole) error("Attempt to drift to the past");
#endif

  /* Drift from the last time the cell was drifted to the current time */
  double dt_drift;
  if (e->policy & engine_policy_cosmology)
    dt_drift =
        cosmology_get_drift_factor(e->cosmology, ti_old_multipole, ti_current);
  else
    dt_drift = (ti_current - ti_old_multipole) * e->time_base;

  /* Drift the multipole */
  if (ti_current > ti_old_multipole) gravity_drift(c->grav.multipole, dt_drift);

  /* Are we not in a leaf ? */
  if (c->split) {
    /* Loop over the progeny and recurse. */
    for (int k = 0; k < 8; k++)
      if (c->progeny[k] != NULL) cell_drift_all_multipoles(c->progeny[k], e);
  }

  /* Update the time of the last drift */
  c->grav.ti_old_multipole = ti_current;
}

/**
 * @brief Drifts the multipole of a cell to the current time.
 *
 * Only drifts the multipole at this level. Multipoles deeper in the
 * tree are not updated.
 *
 * @param c The #cell.
 * @param e The #engine (to get ti_current).
 */
void cell_drift_multipole(struct cell *c, const struct engine *e) {
  const integertime_t ti_old_multipole = c->grav.ti_old_multipole;
  const integertime_t ti_current = e->ti_current;

#ifdef SWIFT_DEBUG_CHECKS
  /* Check that we are actually going to move forward. */
  if (ti_current < ti_old_multipole) error("Attempt to drift to the past");
#endif

  /* Drift from the last time the cell was drifted to the current time */
  double dt_drift;
  if (e->policy & engine_policy_cosmology)
    dt_drift =
        cosmology_get_drift_factor(e->cosmology, ti_old_multipole, ti_current);
  else
    dt_drift = (ti_current - ti_old_multipole) * e->time_base;

  if (ti_current > ti_old_multipole) gravity_drift(c->grav.multipole, dt_drift);

  /* Update the time of the last drift */
  c->grav.ti_old_multipole = ti_current;
}
