/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (schaller@strw.leidenuniv.nl)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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

/* Standard headers. */
#include <stdint.h>

/* This object's header. */
#include "space.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "debug.h"
#include "engine.h"
#include "multipole.h"
#include "star_formation_logger.h"
#include "threadpool.h"

/**
 * @brief Finalise a leaf cell.
 *
 * This function collects the time-step and smoothing length information from
 * the particles in the cell and updates the cell's properties accordingly. This
 * information can then be handed back up the tree.
 *
 * @param s The #space the cell lives in.
 * @param c The leaf #cell to finalise.
 */
static void space_split_finalise_leaf(struct space *s, struct cell *c) {

  /* Unpack cell information. */
  const int count = c->hydro.count;
  const int gcount = c->grav.count;
  const int scount = c->stars.count;
  const int bcount = c->black_holes.count;
  const int sink_count = c->sinks.count;
  struct part *parts = c->hydro.parts;
  struct gpart *gparts = c->grav.parts;
  struct spart *sparts = c->stars.parts;
  struct bpart *bparts = c->black_holes.parts;
  struct xpart *xparts = c->hydro.xparts;
  struct sink *sinks = c->sinks.parts;
  struct engine *e = s->e;
  const integertime_t ti_current = e->ti_current;
  const int with_rt = e->policy & engine_policy_rt;

  /* Initialise the variables we will aggregate. */
  float h_max = 0.f;
  float h_max_active = 0.f;
  float stars_h_max = 0.f;
  float stars_h_max_active = 0.f;
  float black_holes_h_max = 0.f;
  float black_holes_h_max_active = 0.f;
  float sinks_h_max = 0.f;
  float sinks_h_max_active = 0.f;
  integertime_t ti_hydro_end_min = max_nr_timesteps, ti_hydro_beg_max = 0;
  integertime_t ti_rt_end_min = max_nr_timesteps, ti_rt_beg_max = 0;
  integertime_t ti_rt_min_step_size = max_nr_timesteps;
  integertime_t ti_gravity_end_min = max_nr_timesteps, ti_gravity_beg_max = 0;
  integertime_t ti_stars_end_min = max_nr_timesteps, ti_stars_beg_max = 0;
  integertime_t ti_sinks_end_min = max_nr_timesteps, ti_sinks_beg_max = 0;
  integertime_t ti_black_holes_end_min = max_nr_timesteps,
                ti_black_holes_beg_max = 0;

  /* Clear the progeny. */
  bzero(c->progeny, sizeof(struct cell *) * 8);

  /* We are a leaf cell, so we are not split. */
  c->split = 0;

  /* hydro: Get dt_min/dt_max. */
  for (int k = 0; k < count; k++) {
#ifdef SWIFT_DEBUG_CHECKS
    if (parts[k].time_bin == time_bin_not_created)
      error("Extra particle present in space_split()");
    if (parts[k].time_bin == time_bin_inhibited)
      error("Inhibited particle present in space_split()");
#endif

    /* When does this particle's time-step start and end? */
    const timebin_t time_bin = parts[k].time_bin;
    const timebin_t time_bin_rt = parts[k].rt_time_data.time_bin;
    const integertime_t ti_end = get_integer_time_end(ti_current, time_bin);
    const integertime_t ti_beg = get_integer_time_begin(ti_current, time_bin);
    ti_hydro_end_min = min(ti_hydro_end_min, ti_end);
    ti_hydro_beg_max = max(ti_hydro_beg_max, ti_beg);

    if (with_rt) {
      /* Contrary to other physics, RT doesn't have its own particle type.
       * So collect time step data from particles only when we're running
       * with RT. Otherwise, we may find cells which are active or in
       * impossible timezones. Skipping this check results in cells having
       * RT times = max_nr_timesteps or zero, respecively. */
      const integertime_t ti_rt_end =
          get_integer_time_end(ti_current, time_bin_rt);
      const integertime_t ti_rt_beg =
          get_integer_time_begin(ti_current, time_bin_rt);
      const integertime_t ti_rt_step = get_integer_timestep(time_bin_rt);
      ti_rt_end_min = min(ti_rt_end_min, ti_rt_end);
      ti_rt_beg_max = max(ti_rt_beg_max, ti_rt_beg);
      ti_rt_min_step_size = min(ti_rt_min_step_size, ti_rt_step);
    }

    /* Get the maximum smoothing length. */
    h_max = max(h_max, parts[k].h);

    /* Get the maximum smoothing length of active particles. */
    if (part_is_active(&parts[k], e))
      h_max_active = max(h_max_active, parts[k].h);

    /* Set the depth of the particle in the cell. */
    cell_set_part_h_depth(&parts[k], c);

    /* Collect SFR from the particles after rebuilt */
    star_formation_logger_log_inactive_part(&parts[k], &xparts[k],
                                            &c->stars.sfh);
  }

  /* xparts: Reset x_diff */
  for (int k = 0; k < count; k++) {
    xparts[k].x_diff[0] = 0.f;
    xparts[k].x_diff[1] = 0.f;
    xparts[k].x_diff[2] = 0.f;
  }

  /* gparts: Get dt_min/dt_max. */
  for (int k = 0; k < gcount; k++) {
#ifdef SWIFT_DEBUG_CHECKS
    if (gparts[k].time_bin == time_bin_not_created)
      error("Extra g-particle present in space_split()");
    if (gparts[k].time_bin == time_bin_inhibited)
      error("Inhibited g-particle present in space_split()");
#endif

    /* When does this particle's time-step start and end? */
    const timebin_t time_bin = gparts[k].time_bin;
    const integertime_t ti_end = get_integer_time_end(ti_current, time_bin);
    const integertime_t ti_beg = get_integer_time_begin(ti_current, time_bin);
    ti_gravity_end_min = min(ti_gravity_end_min, ti_end);
    ti_gravity_beg_max = max(ti_gravity_beg_max, ti_beg);
  }

  /* sparts: Get dt_min/dt_max */
  for (int k = 0; k < scount; k++) {
#ifdef SWIFT_DEBUG_CHECKS
    if (sparts[k].time_bin == time_bin_not_created)
      error("Extra s-particle present in space_split()");
    if (sparts[k].time_bin == time_bin_inhibited)
      error("Inhibited s-particle present in space_split()");
#endif

    /* When does this particle's time-step start and end? */
    const timebin_t time_bin = sparts[k].time_bin;
    const integertime_t ti_end = get_integer_time_end(ti_current, time_bin);
    const integertime_t ti_beg = get_integer_time_begin(ti_current, time_bin);
    ti_stars_end_min = min(ti_stars_end_min, ti_end);
    ti_stars_beg_max = max(ti_stars_beg_max, ti_beg);

    /* Get the maximum smoothing length. */
    stars_h_max = max(stars_h_max, sparts[k].h);

    /* Get the maximum smoothing length of active particles. */
    if (spart_is_active(&sparts[k], e))
      stars_h_max_active = max(stars_h_max_active, sparts[k].h);

    /* Set the depth of the particle in the cell. */
    cell_set_spart_h_depth(&sparts[k], c);

    /* Reset x_diff */
    sparts[k].x_diff[0] = 0.f;
    sparts[k].x_diff[1] = 0.f;
    sparts[k].x_diff[2] = 0.f;
  }

  /* sinks: Get dt_min/dt_max */
  for (int k = 0; k < sink_count; k++) {
#ifdef SWIFT_DEBUG_CHECKS
    if (sinks[k].time_bin == time_bin_not_created)
      error("Extra sink-particle present in space_split()");
    if (sinks[k].time_bin == time_bin_inhibited)
      error("Inhibited sink-particle present in space_split()");
#endif

    /* When does this particle's time-step start and end? */
    const timebin_t time_bin = sinks[k].time_bin;
    const integertime_t ti_end = get_integer_time_end(ti_current, time_bin);
    const integertime_t ti_beg = get_integer_time_begin(ti_current, time_bin);
    ti_sinks_end_min = min(ti_sinks_end_min, ti_end);
    ti_sinks_beg_max = max(ti_sinks_beg_max, ti_beg);

    /* Get the maximum smoothing length. */
    sinks_h_max = max(sinks_h_max, sinks[k].h);

    /* Get the maximum smoothing length of active particles. */
    if (sink_is_active(&sinks[k], e))
      sinks_h_max_active = max(sinks_h_max_active, sinks[k].h);

    /* Set the depth of the particle in the cell. */
    cell_set_sink_h_depth(&sinks[k], c);

    /* Collect SFR from the particles after rebuilt */
    star_formation_logger_log_inactive_sink(&sinks[k], &c->stars.sfh);

    /* Reset x_diff */
    sinks[k].x_diff[0] = 0.f;
    sinks[k].x_diff[1] = 0.f;
    sinks[k].x_diff[2] = 0.f;
  }

  /* bparts: Get dt_min/dt_max */
  for (int k = 0; k < bcount; k++) {
#ifdef SWIFT_DEBUG_CHECKS
    if (bparts[k].time_bin == time_bin_not_created)
      error("Extra b-particle present in space_split()");
    if (bparts[k].time_bin == time_bin_inhibited)
      error("Inhibited b-particle present in space_split()");
#endif

    /* When does this particle's time-step start and end? */
    const timebin_t time_bin = bparts[k].time_bin;
    const integertime_t ti_end = get_integer_time_end(ti_current, time_bin);
    const integertime_t ti_beg = get_integer_time_begin(ti_current, time_bin);
    ti_black_holes_end_min = min(ti_black_holes_end_min, ti_end);
    ti_black_holes_beg_max = max(ti_black_holes_beg_max, ti_beg);

    /* Get the maximum smoothing length. */
    black_holes_h_max = max(black_holes_h_max, bparts[k].h);

    /* Get the maximum smoothing length of active particles. */
    if (bpart_is_active(&bparts[k], e))
      black_holes_h_max_active = max(black_holes_h_max_active, bparts[k].h);

    /* Set the depth of the particle in the cell. */
    cell_set_bpart_h_depth(&bparts[k], c);

    /* Reset x_diff */
    bparts[k].x_diff[0] = 0.f;
    bparts[k].x_diff[1] = 0.f;
    bparts[k].x_diff[2] = 0.f;
  }

  /* Construct the multipole and the centre of mass*/
  if (s->with_self_gravity) {
    if (gcount > 0) {

      gravity_P2M(c->grav.multipole, c->grav.parts, c->grav.count,
                  e->gravity_properties);

      /* Compute the multipole power */
      gravity_multipole_compute_power(&c->grav.multipole->m_pole);

    } else {

      /* No gparts in that leaf cell */

      /* Set the values to something sensible */
      gravity_multipole_init(&c->grav.multipole->m_pole);
      if (c->nodeID == engine_rank) {
        c->grav.multipole->CoM[0] = c->loc[0] + c->width[0] / 2.;
        c->grav.multipole->CoM[1] = c->loc[1] + c->width[1] / 2.;
        c->grav.multipole->CoM[2] = c->loc[2] + c->width[2] / 2.;
        c->grav.multipole->r_max = 0.;
      }
    }

    /* Store the value at rebuild time */
    c->grav.multipole->r_max_rebuild = c->grav.multipole->r_max;
    c->grav.multipole->CoM_rebuild[0] = c->grav.multipole->CoM[0];
    c->grav.multipole->CoM_rebuild[1] = c->grav.multipole->CoM[1];
    c->grav.multipole->CoM_rebuild[2] = c->grav.multipole->CoM[2];
    c->grav.multipole->dx_max[0] = 0.f;
    c->grav.multipole->dx_max[1] = 0.f;
    c->grav.multipole->dx_max[2] = 0.f;
  }

  /* Set the values for this cell. */
  c->hydro.h_max = h_max;
  c->hydro.h_max_active = h_max_active;
  c->hydro.ti_end_min = ti_hydro_end_min;
  c->hydro.ti_beg_max = ti_hydro_beg_max;
  c->rt.ti_rt_end_min = ti_rt_end_min;
  c->rt.ti_rt_beg_max = ti_rt_beg_max;
  c->rt.ti_rt_min_step_size = ti_rt_min_step_size;
  c->grav.ti_end_min = ti_gravity_end_min;
  c->grav.ti_beg_max = ti_gravity_beg_max;
  c->stars.ti_end_min = ti_stars_end_min;
  c->stars.ti_beg_max = ti_stars_beg_max;
  c->stars.h_max = stars_h_max;
  c->stars.h_max_active = stars_h_max_active;
  c->sinks.ti_end_min = ti_sinks_end_min;
  c->sinks.ti_beg_max = ti_sinks_beg_max;
  c->sinks.h_max = sinks_h_max;
  c->sinks.h_max_active = sinks_h_max_active;
  c->black_holes.ti_end_min = ti_black_holes_end_min;
  c->black_holes.ti_beg_max = ti_black_holes_beg_max;
  c->black_holes.h_max = black_holes_h_max;
  c->black_holes.h_max_active = black_holes_h_max_active;
  c->maxdepth = c->depth;

  /* No runner owns this cell yet. We assign those during scheduling. */
  c->owner = -1;
}

/**
 * @brief Accumulate the cell level particle properties from progeny.
 *
 * This takes the max/min of c's current fields against cp's, so c's relevant
 * fields must already be reset to the min/max for the reduction (0 for the
 * h_max fields, max_nr_timesteps for the *_end_min and ti_rt_min_step_size
 * fields, 0 for the *_beg_max fields) before calling this function.
 *
 * @param c The parent #cell, updated in place.
 * @param cp The child #cell being folded in.
 */
static void space_split_accumulate_props(struct cell *c,
                                         const struct cell *cp) {

  /* Smoothing lengths */
  c->hydro.h_max = max(c->hydro.h_max, cp->hydro.h_max);
  c->hydro.h_max_active = max(c->hydro.h_max_active, cp->hydro.h_max_active);
  c->stars.h_max = max(c->stars.h_max, cp->stars.h_max);
  c->stars.h_max_active = max(c->stars.h_max_active, cp->stars.h_max_active);
  c->black_holes.h_max = max(c->black_holes.h_max, cp->black_holes.h_max);
  c->black_holes.h_max_active =
      max(c->black_holes.h_max_active, cp->black_holes.h_max_active);
  c->sinks.h_max = max(c->sinks.h_max, cp->sinks.h_max);
  c->sinks.h_max_active = max(c->sinks.h_max_active, cp->sinks.h_max_active);

  /* Time-step information */
  c->hydro.ti_end_min = min(c->hydro.ti_end_min, cp->hydro.ti_end_min);
  c->hydro.ti_beg_max = max(c->hydro.ti_beg_max, cp->hydro.ti_beg_max);
  c->rt.ti_rt_end_min = min(c->rt.ti_rt_end_min, cp->rt.ti_rt_end_min);
  c->rt.ti_rt_beg_max = max(c->rt.ti_rt_beg_max, cp->rt.ti_rt_beg_max);
  c->rt.ti_rt_min_step_size =
      min(c->rt.ti_rt_min_step_size, cp->rt.ti_rt_min_step_size);
  c->grav.ti_end_min = min(c->grav.ti_end_min, cp->grav.ti_end_min);
  c->grav.ti_beg_max = max(c->grav.ti_beg_max, cp->grav.ti_beg_max);
  c->stars.ti_end_min = min(c->stars.ti_end_min, cp->stars.ti_end_min);
  c->stars.ti_beg_max = max(c->stars.ti_beg_max, cp->stars.ti_beg_max);
  c->sinks.ti_end_min = min(c->sinks.ti_end_min, cp->sinks.ti_end_min);
  c->sinks.ti_beg_max = max(c->sinks.ti_beg_max, cp->sinks.ti_beg_max);
  c->black_holes.ti_end_min =
      min(c->black_holes.ti_end_min, cp->black_holes.ti_end_min);
  c->black_holes.ti_beg_max =
      max(c->black_holes.ti_beg_max, cp->black_holes.ti_beg_max);

  /* Star formation history */
  star_formation_logger_add(&c->stars.sfh, &cp->stars.sfh);
}

/**
 * @brief Build a cell's multipole from its progeny's multipoles (M2M).
 *
 * This function populates a cell's multipole based on the progeny multipoles.
 * These multipoles are constructed from the bottom up with the leaf cells
 * populated based on the gparts they hold.
 *
 * @param c The parent #cell whose multipole is (re)built from its progeny.
 */
static void space_split_populate_multipole(struct cell *c) {

  /* Reset everything */
  gravity_reset(c->grav.multipole);

  /* Compute CoM and bulk velocity from all progenies */
  double CoM[3] = {0., 0., 0.};
  double vel[3] = {0., 0., 0.};
  float max_delta_vel[3] = {0.f, 0.f, 0.f};
  float min_delta_vel[3] = {0.f, 0.f, 0.f};
  double mass = 0.;

  for (int k = 0; k < 8; ++k) {
    if (c->progeny[k] != NULL) {
      const struct gravity_tensors *m = c->progeny[k]->grav.multipole;

      mass += m->m_pole.M_000;

      CoM[0] += m->CoM[0] * m->m_pole.M_000;
      CoM[1] += m->CoM[1] * m->m_pole.M_000;
      CoM[2] += m->CoM[2] * m->m_pole.M_000;

      vel[0] += m->m_pole.vel[0] * m->m_pole.M_000;
      vel[1] += m->m_pole.vel[1] * m->m_pole.M_000;
      vel[2] += m->m_pole.vel[2] * m->m_pole.M_000;

      max_delta_vel[0] = max(m->m_pole.max_delta_vel[0], max_delta_vel[0]);
      max_delta_vel[1] = max(m->m_pole.max_delta_vel[1], max_delta_vel[1]);
      max_delta_vel[2] = max(m->m_pole.max_delta_vel[2], max_delta_vel[2]);

      min_delta_vel[0] = min(m->m_pole.min_delta_vel[0], min_delta_vel[0]);
      min_delta_vel[1] = min(m->m_pole.min_delta_vel[1], min_delta_vel[1]);
      min_delta_vel[2] = min(m->m_pole.min_delta_vel[2], min_delta_vel[2]);
    }
  }

  /* Final operation on the CoM and bulk velocity */
  const double inv_mass = 1. / mass;
  c->grav.multipole->CoM[0] = CoM[0] * inv_mass;
  c->grav.multipole->CoM[1] = CoM[1] * inv_mass;
  c->grav.multipole->CoM[2] = CoM[2] * inv_mass;
  c->grav.multipole->m_pole.vel[0] = vel[0] * inv_mass;
  c->grav.multipole->m_pole.vel[1] = vel[1] * inv_mass;
  c->grav.multipole->m_pole.vel[2] = vel[2] * inv_mass;

  /* Min max velocity along each axis */
  c->grav.multipole->m_pole.max_delta_vel[0] = max_delta_vel[0];
  c->grav.multipole->m_pole.max_delta_vel[1] = max_delta_vel[1];
  c->grav.multipole->m_pole.max_delta_vel[2] = max_delta_vel[2];
  c->grav.multipole->m_pole.min_delta_vel[0] = min_delta_vel[0];
  c->grav.multipole->m_pole.min_delta_vel[1] = min_delta_vel[1];
  c->grav.multipole->m_pole.min_delta_vel[2] = min_delta_vel[2];

  /* Now shift progeny multipoles and add them up */
  struct multipole temp;
  double r_max = 0.;
  for (int k = 0; k < 8; ++k) {
    if (c->progeny[k] != NULL) {
      const struct cell *cp = c->progeny[k];
      const struct multipole *m = &cp->grav.multipole->m_pole;

      /* Contribution to multipole */
      gravity_M2M(&temp, m, c->grav.multipole->CoM, cp->grav.multipole->CoM);
      gravity_multipole_add(&c->grav.multipole->m_pole, &temp);

      /* Upper limit of max CoM<->gpart distance */
      const double dx = c->grav.multipole->CoM[0] - cp->grav.multipole->CoM[0];
      const double dy = c->grav.multipole->CoM[1] - cp->grav.multipole->CoM[1];
      const double dz = c->grav.multipole->CoM[2] - cp->grav.multipole->CoM[2];
      const double r2 = dx * dx + dy * dy + dz * dz;
      r_max = max(r_max, cp->grav.multipole->r_max + sqrt(r2));
    }
  }

  /* Alternative upper limit of max CoM<->gpart distance */
  const double dx = c->grav.multipole->CoM[0] > c->loc[0] + c->width[0] / 2.
                        ? c->grav.multipole->CoM[0] - c->loc[0]
                        : c->loc[0] + c->width[0] - c->grav.multipole->CoM[0];
  const double dy = c->grav.multipole->CoM[1] > c->loc[1] + c->width[1] / 2.
                        ? c->grav.multipole->CoM[1] - c->loc[1]
                        : c->loc[1] + c->width[1] - c->grav.multipole->CoM[1];
  const double dz = c->grav.multipole->CoM[2] > c->loc[2] + c->width[2] / 2.
                        ? c->grav.multipole->CoM[2] - c->loc[2]
                        : c->loc[2] + c->width[2] - c->grav.multipole->CoM[2];

  /* Take minimum of both limits */
  c->grav.multipole->r_max = min(r_max, sqrt(dx * dx + dy * dy + dz * dz));

  /* Store the value at rebuild time */
  c->grav.multipole->r_max_rebuild = c->grav.multipole->r_max;
  c->grav.multipole->CoM_rebuild[0] = c->grav.multipole->CoM[0];
  c->grav.multipole->CoM_rebuild[1] = c->grav.multipole->CoM[1];
  c->grav.multipole->CoM_rebuild[2] = c->grav.multipole->CoM[2];
  c->grav.multipole->dx_max[0] = 0.f;
  c->grav.multipole->dx_max[1] = 0.f;
  c->grav.multipole->dx_max[2] = 0.f;

  /* Compute the multipole power */
  gravity_multipole_compute_power(&c->grav.multipole->m_pole);
}

/**
 * @brief Fill a top-level cell's slice of the space-wide sorting buffers.
 *
 * The buffers are allocated once, up front, for the complete local
 * particle arrays by #space_split() -- see #space_split_mapper() -- rather
 * than per top-level cell. This only fills a single cell's slice of those
 * buffers from its own particles; buff, gbuff, sbuff, bbuff and sink_buff
 * must already point at the start of that cell's slice (i.e. offset by the
 * same particle-array offset used elsewhere, e.g. in #cell_split()), and
 * the *_offset arguments must be that same offset, so that each entry's
 * part_ind can be set to the particle's index in the space-wide array
 * (e.g. s->parts) it was filled from. This is what lets #cell_split()
 * reorder the buffers alone, further down the tree, without also having to
 * move the particles themselves.
 *
 * @param c The #cell whose particles populate the buffers.
 * @param buff This cell's slice of the space-wide hydro buffer.
 * @param gbuff This cell's slice of the space-wide gravity buffer.
 * @param sbuff This cell's slice of the space-wide star buffer.
 * @param bbuff This cell's slice of the space-wide black hole buffer.
 * @param sink_buff This cell's slice of the space-wide sink buffer.
 * @param parts_offset c->hydro.parts - s->parts.
 * @param gparts_offset c->grav.parts - s->gparts.
 * @param sparts_offset c->stars.parts - s->sparts.
 * @param bparts_offset c->black_holes.parts - s->bparts.
 * @param sinks_offset c->sinks.parts - s->sinks.
 */
static void space_split_fill_buffers(
    struct cell *c, struct cell_buff *restrict buff,
    struct cell_buff *restrict gbuff, struct cell_buff *restrict sbuff,
    struct cell_buff *restrict bbuff, struct cell_buff *restrict sink_buff,
    const ptrdiff_t parts_offset, const ptrdiff_t gparts_offset,
    const ptrdiff_t sparts_offset, const ptrdiff_t bparts_offset,
    const ptrdiff_t sinks_offset) {

  /* Unpack counts and particle arrays. */
  const int count = c->hydro.count;
  const int gcount = c->grav.count;
  const int scount = c->stars.count;
  const int bcount = c->black_holes.count;
  const int sink_count = c->sinks.count;
  struct part *parts = c->hydro.parts;
  struct gpart *gparts = c->grav.parts;
  struct spart *sparts = c->stars.parts;
  struct bpart *bparts = c->black_holes.parts;
  struct sink *sinks = c->sinks.parts;

  /* Fill the temporary buffer for hydro parts. */
  for (int k = 0; k < count; k++) {
#ifdef SWIFT_DEBUG_CHECKS
    if (parts[k].time_bin == time_bin_inhibited)
      error("Inhibited particle present in space_split()");
    if (parts[k].time_bin == time_bin_not_created)
      error("Extra particle present in space_split()");
#endif
    buff[k].x[0] = parts[k].x[0];
    buff[k].x[1] = parts[k].x[1];
    buff[k].x[2] = parts[k].x[2];
    buff[k].part_ind = parts_offset + k;
  }

  /* Fill the temporary buffer for gravity parts. */
  for (int k = 0; k < gcount; k++) {
#ifdef SWIFT_DEBUG_CHECKS
    if (gparts[k].time_bin == time_bin_inhibited)
      error("Inhibited particle present in space_split()");
    if (gparts[k].time_bin == time_bin_not_created)
      error("Extra particle present in space_split()");
#endif
    gbuff[k].x[0] = gparts[k].x[0];
    gbuff[k].x[1] = gparts[k].x[1];
    gbuff[k].x[2] = gparts[k].x[2];
    gbuff[k].part_ind = gparts_offset + k;
  }

  /* Fill the temporary buffer for star parts. */
  for (int k = 0; k < scount; k++) {
#ifdef SWIFT_DEBUG_CHECKS
    if (sparts[k].time_bin == time_bin_inhibited)
      error("Inhibited particle present in space_split()");
    if (sparts[k].time_bin == time_bin_not_created)
      error("Extra particle present in space_split()");
#endif
    sbuff[k].x[0] = sparts[k].x[0];
    sbuff[k].x[1] = sparts[k].x[1];
    sbuff[k].x[2] = sparts[k].x[2];
    sbuff[k].part_ind = sparts_offset + k;
  }

  /* Fill the temporary buffer for black hole parts. */
  for (int k = 0; k < bcount; k++) {
#ifdef SWIFT_DEBUG_CHECKS
    if (bparts[k].time_bin == time_bin_inhibited)
      error("Inhibited particle present in space_split()");
    if (bparts[k].time_bin == time_bin_not_created)
      error("Extra particle present in space_split()");
#endif
    bbuff[k].x[0] = bparts[k].x[0];
    bbuff[k].x[1] = bparts[k].x[1];
    bbuff[k].x[2] = bparts[k].x[2];
    bbuff[k].part_ind = bparts_offset + k;
  }

  /* Fill the temporary buffer for sink parts. */
  for (int k = 0; k < sink_count; k++) {
#ifdef SWIFT_DEBUG_CHECKS
    if (sinks[k].time_bin == time_bin_inhibited)
      error("Inhibited particle present in space_split()");
    if (sinks[k].time_bin == time_bin_not_created)
      error("Extra particle present in space_split()");
#endif
    sink_buff[k].x[0] = sinks[k].x[0];
    sink_buff[k].x[1] = sinks[k].x[1];
    sink_buff[k].x[2] = sinks[k].x[2];
    sink_buff[k].part_ind = sinks_offset + k;
  }
}

/**
 * @brief A thread's reused, per-top-level-cell particle scratch buffers.
 *
 * Sized to the largest top-level cell this thread has processed so far
 * (grown, never shrunk, across the whole run -- see
 * #space_split_ensure_particle_scratch()), not to the space's complete
 * local particle arrays: #space_split_move_leaf() gathers each leaf's
 * particles into these, and once a whole top-level cell's subtree has been
 * gathered, #space_split_mapper() copies each family's slice straight back
 * into the real array with a single bulk memcpy. That keeps the one real
 * per-particle data move this design allows for (see
 * #space_split_move_leaf()) sequential on both ends, without ever needing
 * an allocation sized to the whole space.
 */
struct space_split_particle_scratch {

  /*! Scratch #part buffer. */
  struct part *parts;

  /*! Scratch #xpart buffer; always sized alongside #parts. */
  struct xpart *xparts;

  /*! Allocated capacity of #parts and #xparts. */
  int parts_capacity;

  /*! Scratch #gpart buffer. */
  struct gpart *gparts;

  /*! Allocated capacity of #gparts. */
  int gparts_capacity;

  /*! Scratch #spart buffer. */
  struct spart *sparts;

  /*! Allocated capacity of #sparts. */
  int sparts_capacity;

  /*! Scratch #bpart buffer. */
  struct bpart *bparts;

  /*! Allocated capacity of #bparts. */
  int bparts_capacity;

  /*! Scratch #sink buffer. */
  struct sink *sinks;

  /*! Allocated capacity of #sinks. */
  int sinks_capacity;
};

/**
 * @brief Grow a thread's particle scratch buffers if they are not already
 *        big enough for a top-level cell with the given per-family counts.
 *
 * Each family is grown independently, since a top-level cell's counts for
 * different families are generally not the same. Growth never needs to
 * preserve a buffer's previous contents -- every entry is fully overwritten
 * by the gather that follows -- so this simply frees the too-small buffer
 * and allocates a bigger one.
 *
 * @param scratch The thread's #space_split_particle_scratch.
 * @param count Hydro particle count required.
 * @param gcount Gravity particle count required.
 * @param scount Star particle count required.
 * @param bcount Black hole particle count required.
 * @param sink_count Sink particle count required.
 */
static void space_split_ensure_particle_scratch(
    struct space_split_particle_scratch *scratch, const int count,
    const int gcount, const int scount, const int bcount,
    const int sink_count) {

  if (count > scratch->parts_capacity) {
    if (scratch->parts != NULL) swift_free("particle_scratch", scratch->parts);
    if (scratch->xparts != NULL)
      swift_free("particle_scratch", scratch->xparts);
    if (swift_memalign("particle_scratch", (void **)&scratch->parts,
                       part_align, sizeof(struct part) * count) != 0)
      error("Failed to allocate particle scratch parts buffer.");
    if (swift_memalign("particle_scratch", (void **)&scratch->xparts,
                       xpart_align, sizeof(struct xpart) * count) != 0)
      error("Failed to allocate particle scratch xparts buffer.");
    scratch->parts_capacity = count;
  }
  if (gcount > scratch->gparts_capacity) {
    if (scratch->gparts != NULL)
      swift_free("particle_scratch", scratch->gparts);
    if (swift_memalign("particle_scratch", (void **)&scratch->gparts,
                       gpart_align, sizeof(struct gpart) * gcount) != 0)
      error("Failed to allocate particle scratch gparts buffer.");
    scratch->gparts_capacity = gcount;
  }
  if (scount > scratch->sparts_capacity) {
    if (scratch->sparts != NULL)
      swift_free("particle_scratch", scratch->sparts);
    if (swift_memalign("particle_scratch", (void **)&scratch->sparts,
                       spart_align, sizeof(struct spart) * scount) != 0)
      error("Failed to allocate particle scratch sparts buffer.");
    scratch->sparts_capacity = scount;
  }
  if (bcount > scratch->bparts_capacity) {
    if (scratch->bparts != NULL)
      swift_free("particle_scratch", scratch->bparts);
    if (swift_memalign("particle_scratch", (void **)&scratch->bparts,
                       bpart_align, sizeof(struct bpart) * bcount) != 0)
      error("Failed to allocate particle scratch bparts buffer.");
    scratch->bparts_capacity = bcount;
  }
  if (sink_count > scratch->sinks_capacity) {
    if (scratch->sinks != NULL)
      swift_free("particle_scratch", scratch->sinks);
    if (swift_memalign("particle_scratch", (void **)&scratch->sinks,
                       sink_align, sizeof(struct sink) * sink_count) != 0)
      error("Failed to allocate particle scratch sinks buffer.");
    scratch->sinks_capacity = sink_count;
  }
}

/**
 * @brief Free a thread's particle scratch buffers.
 *
 * NULL-guarded and NULL-setting, so it is idempotent: it is called once per
 * worker thread from #space_split_aggregate_mapper() (freeing the buffers
 * in parallel, on the thread and NUMA node that allocated and touched them,
 * rather than serially on the main thread), and may be called again for the
 * same thread if that thread handles more than one aggregation chunk, and
 * finally by #space_split() as a serial straggler sweep for any thread that
 * was never handed an aggregation chunk. All but the first call for a given
 * thread are no-ops.
 *
 * @param scratch The thread's #space_split_particle_scratch.
 */
static void space_split_free_particle_scratch(
    struct space_split_particle_scratch *scratch) {

  if (scratch->parts != NULL) {
    swift_free("particle_scratch", scratch->parts);
    scratch->parts = NULL;
  }
  if (scratch->xparts != NULL) {
    swift_free("particle_scratch", scratch->xparts);
    scratch->xparts = NULL;
  }
  if (scratch->gparts != NULL) {
    swift_free("particle_scratch", scratch->gparts);
    scratch->gparts = NULL;
  }
  if (scratch->sparts != NULL) {
    swift_free("particle_scratch", scratch->sparts);
    scratch->sparts = NULL;
  }
  if (scratch->bparts != NULL) {
    swift_free("particle_scratch", scratch->bparts);
    scratch->bparts = NULL;
  }
  if (scratch->sinks != NULL) {
    swift_free("particle_scratch", scratch->sinks);
    scratch->sinks = NULL;
  }
}

/**
 * @brief Gather a leaf cell's particles into its slice of the thread's
 *        particle scratch buffers.
 *
 * Called from #space_split_recursive() the moment a cell is found not to
 * need any further splitting, while its buffer slice (buff/sbuff/bbuff/
 * gbuff/sink_buff, already fully sorted for this cell by the parent's
 * #cell_split() calls) is still hot. For each family, and each entry k in
 * the leaf's slice, this reads the source particle identified by part_ind
 * and writes it into the matching slot of this leaf's slice of the scratch
 * buffers -- @p scratch_parts etc. have already been advanced by the
 * caller to point at this leaf's own position within the top-level cell's
 * scratch buffers, exactly like @p buff etc. have for the sort buffers, so
 * no further offset is needed to address them.
 *
 * A particle already in its final slot (source index == final index) is
 * *not* gathered: it stays put in the real array and its scratch slot is
 * left untouched, to be skipped on write-back (see #space_split_mapper()).
 * This is the same work master avoids for particles that have not changed
 * cell. Its partner gpart's back-index is still fixed up, though, since
 * that gpart may itself be moving.
 *
 * #gpart linkage is fixed up as each particle lands, exactly as it would
 * be in a single global pass, but scoped to this one leaf: a #gpart and
 * its #part/#spart/#bpart/#sink partner always share identical
 * coordinates, so #cell_split()'s consistent (`>=`) pivot comparison
 * guarantees they always compute the same bucket at every level and
 * therefore always end up in the same leaf. That is what makes it safe to
 * do all of this leaf's relinking self-contained, right here: the
 * non-gravity families are gathered first, each fixing up its partner
 * #gpart's (not yet gathered, so still safe to reach at its old, real
 * location) id_or_neg_offset; #gpart is gathered last, using that
 * now-correct id_or_neg_offset to fix up the `.gpart` back-pointer to its
 * own real, final address. That back-pointer must land wherever write-back
 * will read the partner baryon from: its scratch slice if the baryon is
 * moving, and -- since a stationary baryon is skipped on write-back -- the
 * real array too if it is not. Writing the real array directly here is
 * safe because a stationary particle is never another leaf's gather source
 * and so is never read concurrently.
 *
 * @param s The #space.
 * @param c The leaf #cell whose particles are gathered.
 * @param scratch_parts This leaf's slice of the thread's scratch #part
 *        buffer.
 * @param scratch_xparts This leaf's slice of the thread's scratch #xpart
 *        buffer.
 * @param scratch_sparts This leaf's slice of the thread's scratch #spart
 *        buffer.
 * @param scratch_bparts This leaf's slice of the thread's scratch #bpart
 *        buffer.
 * @param scratch_sinks This leaf's slice of the thread's scratch #sink
 *        buffer.
 * @param scratch_gparts This leaf's slice of the thread's scratch #gpart
 *        buffer.
 * @param buff This leaf's slice of the space-wide hydro buffer.
 * @param sbuff This leaf's slice of the space-wide star buffer.
 * @param bbuff This leaf's slice of the space-wide black hole buffer.
 * @param gbuff This leaf's slice of the space-wide gravity buffer.
 * @param sink_buff This leaf's slice of the space-wide sink buffer.
 */
static void space_split_move_leaf(
    struct space *s, struct cell *c, struct part *restrict scratch_parts,
    struct xpart *restrict scratch_xparts,
    struct spart *restrict scratch_sparts,
    struct bpart *restrict scratch_bparts,
    struct sink *restrict scratch_sinks,
    struct gpart *restrict scratch_gparts, struct cell_buff *restrict buff,
    struct cell_buff *restrict sbuff, struct cell_buff *restrict bbuff,
    struct cell_buff *restrict gbuff, struct cell_buff *restrict sink_buff) {

  const int count = c->hydro.count;
  const int scount = c->stars.count;
  const int bcount = c->black_holes.count;
  const int sink_count = c->sinks.count;
  const int gcount = c->grav.count;

  /* Hydro: gather parts and xparts, fixing up the id_or_neg_offset of the
   * gpart each is linked to (not yet gathered, so this is safe). A particle
   * already in its final slot (j == i) is left untouched in the real array
   * -- its scratch slot stays garbage and is skipped on write-back -- but
   * its partner gpart's back-index is still fixed up unconditionally, since
   * that gpart may itself be moving. */
  if (count > 0) {
    const ptrdiff_t offset = c->hydro.parts - s->parts;
    for (int k = 0; k < count; k++) {
      const size_t j = buff[k].part_ind;
      const size_t i = offset + k;
      struct part *src = &s->parts[j];
      if (j != i) {
        scratch_parts[k] = *src;
        scratch_xparts[k] = s->xparts[j];
      }
      if (src->gpart != NULL) src->gpart->id_or_neg_offset = -(ptrdiff_t)i;
    }
  }

  /* Stars. */
  if (scount > 0) {
    const ptrdiff_t offset = c->stars.parts - s->sparts;
    for (int k = 0; k < scount; k++) {
      const size_t j = sbuff[k].part_ind;
      const size_t i = offset + k;
      struct spart *src = &s->sparts[j];
      if (j != i) scratch_sparts[k] = *src;
      if (src->gpart != NULL) src->gpart->id_or_neg_offset = -(ptrdiff_t)i;
    }
  }

  /* Black holes. */
  if (bcount > 0) {
    const ptrdiff_t offset = c->black_holes.parts - s->bparts;
    for (int k = 0; k < bcount; k++) {
      const size_t j = bbuff[k].part_ind;
      const size_t i = offset + k;
      struct bpart *src = &s->bparts[j];
      if (j != i) scratch_bparts[k] = *src;
      if (src->gpart != NULL) src->gpart->id_or_neg_offset = -(ptrdiff_t)i;
    }
  }

  /* Sinks. */
  if (sink_count > 0) {
    const ptrdiff_t offset = c->sinks.parts - s->sinks;
    for (int k = 0; k < sink_count; k++) {
      const size_t j = sink_buff[k].part_ind;
      const size_t i = offset + k;
      struct sink *src = &s->sinks[j];
      if (j != i) scratch_sinks[k] = *src;
      if (src->gpart != NULL) src->gpart->id_or_neg_offset = -(ptrdiff_t)i;
    }
  }

  /* Gravity -- last, using the id_or_neg_offset values the four families
   * above just fixed up, to fix up their `.gpart` back-pointers in turn. */
  if (gcount > 0) {
    const ptrdiff_t hydro_offset = c->hydro.parts - s->parts;
    const ptrdiff_t stars_offset = c->stars.parts - s->sparts;
    const ptrdiff_t bh_offset = c->black_holes.parts - s->bparts;
    const ptrdiff_t sink_offset = c->sinks.parts - s->sinks;
    const ptrdiff_t offset = c->grav.parts - s->gparts;
    for (int k = 0; k < gcount; k++) {
      const size_t j = gbuff[k].part_ind;
      const size_t i = offset + k;
      struct gpart *src = &s->gparts[j];
      if (j != i) scratch_gparts[k] = *src;

      /* The real, final address this gpart will occupy once this
       * top-level cell's scratch data is copied back -- valid to store
       * now even though that memory does not hold this data yet.
       *
       * The partner baryon's `.gpart` back-pointer must land wherever
       * write-back will read that baryon from: its scratch slot if the
       * baryon is moving; and, if the baryon is stationary (and hence
       * skipped on write-back), directly in the real array at its final
       * index too. A stationary particle is never another leaf's gather
       * source, so this direct real-array write cannot race any gather. */
      struct gpart *final_gp = &s->gparts[i];
      const ptrdiff_t partner = -src->id_or_neg_offset;
      if (src->type == swift_type_gas) {
        const ptrdiff_t pl = partner - hydro_offset;
        scratch_parts[pl].gpart = final_gp;
        if (buff[pl].part_ind == (size_t)partner)
          s->parts[partner].gpart = final_gp;
      } else if (src->type == swift_type_stars) {
        const ptrdiff_t pl = partner - stars_offset;
        scratch_sparts[pl].gpart = final_gp;
        if (sbuff[pl].part_ind == (size_t)partner)
          s->sparts[partner].gpart = final_gp;
      } else if (src->type == swift_type_black_hole) {
        const ptrdiff_t pl = partner - bh_offset;
        scratch_bparts[pl].gpart = final_gp;
        if (bbuff[pl].part_ind == (size_t)partner)
          s->bparts[partner].gpart = final_gp;
      } else if (src->type == swift_type_sink) {
        const ptrdiff_t pl = partner - sink_offset;
        scratch_sinks[pl].gpart = final_gp;
        if (sink_buff[pl].part_ind == (size_t)partner)
          s->sinks[partner].gpart = final_gp;
      }
    }
  }
}

/**
 * @brief Recursively split a cell.
 *
 * This builds the cell hierarchy and sorts particles into the progeny; it
 * does not compute any per-cell statistics (h_max, time-step bounds, star
 * formation history, multipoles, ...). Those are handled for every cell --
 * leaf and split alike -- by the separate aggregation pass in
 * #space_split_aggregate_recursive(), once the whole hierarchy below a
 * top-level cell has been built here.
 *
 * Leaves are different: the moment a cell is found not to need any further
 * splitting, its particles are gathered into this thread's particle
 * scratch buffers right here, via #space_split_move_leaf(), while its
 * buffer slice is still hot -- see that function's documentation for why
 * this is both correct and the point of the exercise. Once the whole
 * top-level cell's subtree has been gathered this way,
 * #space_split_mapper() copies the scratch buffers back into the real
 * arrays with a single bulk memcpy per family.
 *
 * @param s The #space in which the cell lives.
 * @param c The #cell to split recursively.
 * @param scratch_parts This cell's slice of the thread's scratch #part
 *        buffer, advanced in step with @p buff.
 * @param scratch_xparts This cell's slice of the thread's scratch #xpart
 *        buffer.
 * @param scratch_sparts This cell's slice of the thread's scratch #spart
 *        buffer, advanced in step with @p sbuff.
 * @param scratch_bparts This cell's slice of the thread's scratch #bpart
 *        buffer, advanced in step with @p bbuff.
 * @param scratch_sinks This cell's slice of the thread's scratch #sink
 *        buffer, advanced in step with @p sink_buff.
 * @param scratch_gparts This cell's slice of the thread's scratch #gpart
 *        buffer, advanced in step with @p gbuff.
 * @param ind Scratch space for #cell_split(), reused unchanged at every
 *        level of this top-level cell's recursion (sized by the caller to
 *        the top-level cell's own count, the largest any descendant will
 *        need).
 * @param buff This cell's slice of the space-wide hydro sorting buffer,
 *        already filled (see #space_split_mapper()).
 * @param sbuff This cell's slice of the space-wide star sorting buffer,
 *        already filled.
 * @param bbuff This cell's slice of the space-wide black hole sorting
 *        buffer, already filled.
 * @param gbuff This cell's slice of the space-wide gravity sorting buffer,
 *        already filled.
 * @param sink_buff This cell's slice of the space-wide sink sorting
 *        buffer, already filled.
 */
static void space_split_recursive(
    struct space *s, struct cell *c, struct part *restrict scratch_parts,
    struct xpart *restrict scratch_xparts,
    struct spart *restrict scratch_sparts,
    struct bpart *restrict scratch_bparts,
    struct sink *restrict scratch_sinks,
    struct gpart *restrict scratch_gparts, int *restrict ind,
    struct cell_buff *restrict buff, struct cell_buff *restrict sbuff,
    struct cell_buff *restrict bbuff, struct cell_buff *restrict gbuff,
    struct cell_buff *restrict sink_buff, const short int tpid) {

  /* Unpack cell information. */
  const int count = c->hydro.count;
  const int gcount = c->grav.count;
  const int scount = c->stars.count;
  const int with_self_gravity = s->with_self_gravity;
  const int depth = c->depth;

  /* Set the top level cell tpid. Doing it here ensures top level cells
   * have the same tpid as their progeny. */
  if (depth == 0) c->tpid = tpid;

  /* If the depth is too large, we have a problem and should stop. */
  if (depth > space_cell_maxdepth) {
    error(
        "Exceeded maximum depth (%d) when splitting cells, aborting. This is "
        "most likely due to having too many particles at the exact same "
        "position, making the construction of a tree impossible.",
        space_cell_maxdepth);
  }

  /* Split or let it be? */
  if ((with_self_gravity && gcount > space_splitsize) ||
      (!with_self_gravity &&
       (count > space_splitsize || scount > space_splitsize))) {

    /* No longer just a leaf. */
    c->split = 1;

    /* Create the cell's progeny. */
    space_getcells(s, 8, c->progeny, tpid);
    for (int k = 0; k < 8; k++) {
      struct cell *cp = c->progeny[k];
      cp->hydro.count = 0;
      cp->grav.count = 0;
      cp->stars.count = 0;
      cp->sinks.count = 0;
      cp->black_holes.count = 0;
      cp->hydro.count_total = 0;
      cp->grav.count_total = 0;
      cp->sinks.count_total = 0;
      cp->stars.count_total = 0;
      cp->black_holes.count_total = 0;
      cp->hydro.ti_old_part = c->hydro.ti_old_part;
      cp->grav.ti_old_part = c->grav.ti_old_part;
      cp->grav.ti_old_multipole = c->grav.ti_old_multipole;
      cp->stars.ti_old_part = c->stars.ti_old_part;
      cp->sinks.ti_old_part = c->sinks.ti_old_part;
      cp->black_holes.ti_old_part = c->black_holes.ti_old_part;
      cp->loc[0] = c->loc[0];
      cp->loc[1] = c->loc[1];
      cp->loc[2] = c->loc[2];
      cp->width[0] = c->width[0] / 2;
      cp->width[1] = c->width[1] / 2;
      cp->width[2] = c->width[2] / 2;
      cp->dmin = c->dmin / 2;
      cp->h_min_allowed = cp->dmin * 0.5 * (1. / kernel_gamma);
      cp->h_max_allowed = cp->dmin * (1. / kernel_gamma);
      if (k & 4) cp->loc[0] += cp->width[0];
      if (k & 2) cp->loc[1] += cp->width[1];
      if (k & 1) cp->loc[2] += cp->width[2];
      cp->depth = c->depth + 1;
      cp->split = 0;
      cp->hydro.h_max = 0.f;
      cp->hydro.h_max_active = 0.f;
      cp->hydro.dx_max_part = 0.f;
      cp->hydro.dx_max_sort = 0.f;
      cp->stars.h_max = 0.f;
      cp->stars.h_max_active = 0.f;
      cp->stars.dx_max_part = 0.f;
      cp->stars.dx_max_sort = 0.f;
      cp->sinks.h_max = 0.f;
      cp->sinks.h_max_active = 0.f;
      cp->sinks.dx_max_part = 0.f;
      cp->black_holes.h_max = 0.f;
      cp->black_holes.h_max_active = 0.f;
      cp->black_holes.dx_max_part = 0.f;
      cp->nodeID = c->nodeID;
      cp->parent = c;
      cp->top = c->top;
      cp->super = NULL;
      cp->hydro.super = NULL;
      cp->grav.super = NULL;
      cp->flags = 0;
      star_formation_logger_init(&cp->stars.sfh);
#ifdef WITH_MPI
      cp->mpi.tag = -1;
#endif  // WITH_MPI
#if defined(SWIFT_DEBUG_CHECKS) || defined(SWIFT_CELL_GRAPH)
      cell_assign_cell_index(cp, c);
#endif
    }

    /* Split the cell's particle data. */
    cell_split(c, ind, buff, sbuff, bbuff, gbuff, sink_buff);

    /* Buffers for the progenitors */
    struct cell_buff *progeny_buff = buff, *progeny_gbuff = gbuff,
                     *progeny_sbuff = sbuff, *progeny_bbuff = bbuff,
                     *progeny_sink_buff = sink_buff;
    struct part *progeny_scratch_parts = scratch_parts;
    struct xpart *progeny_scratch_xparts = scratch_xparts;
    struct spart *progeny_scratch_sparts = scratch_sparts;
    struct bpart *progeny_scratch_bparts = scratch_bparts;
    struct sink *progeny_scratch_sinks = scratch_sinks;
    struct gpart *progeny_scratch_gparts = scratch_gparts;

    for (int k = 0; k < 8; k++) {

      /* Get the progenitor */
      struct cell *cp = c->progeny[k];

      /* Remove any progeny with zero particles. */
      if (cp->hydro.count == 0 && cp->grav.count == 0 && cp->stars.count == 0 &&
          cp->black_holes.count == 0 && cp->sinks.count == 0) {

        space_recycle(s, cp);
        c->progeny[k] = NULL;

      } else {

        /* Recurse */
        space_split_recursive(s, cp, progeny_scratch_parts,
                              progeny_scratch_xparts, progeny_scratch_sparts,
                              progeny_scratch_bparts, progeny_scratch_sinks,
                              progeny_scratch_gparts, ind, progeny_buff,
                              progeny_sbuff, progeny_bbuff, progeny_gbuff,
                              progeny_sink_buff, tpid);

        /* Update the pointers in the buffers */
        progeny_buff += cp->hydro.count;
        progeny_gbuff += cp->grav.count;
        progeny_sbuff += cp->stars.count;
        progeny_bbuff += cp->black_holes.count;
        progeny_sink_buff += cp->sinks.count;
        progeny_scratch_parts += cp->hydro.count;
        progeny_scratch_xparts += cp->hydro.count;
        progeny_scratch_gparts += cp->grav.count;
        progeny_scratch_sparts += cp->stars.count;
        progeny_scratch_bparts += cp->black_holes.count;
        progeny_scratch_sinks += cp->sinks.count;
      }
    }

  } /* Split or let it be? */

  /* Otherwise, this cell remains a leaf: mark it as such, and gather its
   * particles into scratch while its buffer slice is still hot. Statistics
   * and multipole construction are still deferred to the aggregation pass
   * in #space_split_aggregate_recursive(), since they are cheap sequential
   * reads once the particles are in final order, and gain nothing from
   * happening here too. */
  else {
    bzero(c->progeny, sizeof(struct cell *) * 8);
    c->split = 0;
    space_split_move_leaf(s, c, scratch_parts, scratch_xparts, scratch_sparts,
                          scratch_bparts, scratch_sinks, scratch_gparts, buff,
                          sbuff, bbuff, gbuff, sink_buff);
  }
}

/**
 * @brief Recursively aggregate progeny properties onto parent cells.
 *
 * Once the split pass in #space_split_recursive() has built the cell
 * hierarchy and sorted particles into it, gathered every leaf's particles
 * into scratch via #space_split_move_leaf(), and #space_split_mapper() has
 * copied each top-level cell's scratch data back into the real arrays,
 * this function derives every cell's statistics in a single depth-first,
 * post-order traversal. Leaf cells (#cell.split unset) are finalised
 * directly from their own (already moved) particles via
 * #space_split_finalise_leaf(). Split cells are built up from their
 * already-aggregated progeny: smoothing lengths, time-step bounds and star
 * formation history are folded in via #space_split_accumulate_props(),
 * and, if running with self-gravity, the multipole is built from the
 * bottom up (M2M) via #space_split_populate_multipole().
 *
 * No pointer rebasing is needed here: every family is copied back into the
 * same array, at the same addresses, it always occupied, so every cell
 * pointer computed during the split pass is still valid.
 *
 * @param s The #space being aggregated.
 * @param c The cell to aggregate.
 */
static void space_split_aggregate_recursive(struct space *s, struct cell *c) {

  /* Leaves are finalised directly from their own particles. */
  if (!c->split) {
    space_split_finalise_leaf(s, c);
    return;
  }

  /* Aggregate progeny first. */
  for (int k = 0; k < 8; k++) {
    if (c->progeny[k] != NULL)
      space_split_aggregate_recursive(s, c->progeny[k]);
  }

  /* Reset to the identity of the max/min reduction that
   * space_split_accumulate_props() performs below, one surviving child at a
   * time. */
  c->hydro.h_max = 0.f;
  c->hydro.h_max_active = 0.f;
  c->stars.h_max = 0.f;
  c->stars.h_max_active = 0.f;
  c->black_holes.h_max = 0.f;
  c->black_holes.h_max_active = 0.f;
  c->sinks.h_max = 0.f;
  c->sinks.h_max_active = 0.f;
  c->hydro.ti_end_min = max_nr_timesteps;
  c->hydro.ti_beg_max = 0;
  c->rt.ti_rt_end_min = max_nr_timesteps;
  c->rt.ti_rt_beg_max = 0;
  c->rt.ti_rt_min_step_size = max_nr_timesteps;
  c->grav.ti_end_min = max_nr_timesteps;
  c->grav.ti_beg_max = 0;
  c->stars.ti_end_min = max_nr_timesteps;
  c->stars.ti_beg_max = 0;
  c->sinks.ti_end_min = max_nr_timesteps;
  c->sinks.ti_beg_max = 0;
  c->black_holes.ti_end_min = max_nr_timesteps;
  c->black_holes.ti_beg_max = 0;
  int maxdepth = 0;

  /* Loop over progeny and accumulate their properties into this cell. */
  for (int k = 0; k < 8; k++) {
    const struct cell *cp = c->progeny[k];
    if (cp == NULL) continue;

    /* Update the cell-wide properties. */
    space_split_accumulate_props(c, cp);

    /* Update the maximum depth. */
    maxdepth = max(maxdepth, cp->maxdepth);
  }

  /* Deal with the multipole */
  if (s->with_self_gravity) space_split_populate_multipole(c);

  /* The per-family summary fields were already written directly into c by
   * space_split_accumulate_props() as each child was folded in. */
  c->maxdepth = maxdepth;

  /* No runner owns this cell yet. We assign those during scheduling. */
  c->owner = -1;
}

/**
 * @brief A thread's reused #cell_split() scratch buffer for bucket ids.
 *
 * See #cell_split()'s documentation for why this is not part of
 * #cell_buff. Grown, never shrunk, as this thread encounters top-level
 * cells with larger counts; reused unchanged for every level of every top
 * cell's recursion, since a descendant's count can only ever be smaller
 * than its top-level cell's own.
 */
struct space_split_ind_scratch {

  /*! The scratch buffer itself. */
  int *ind;

  /*! Number of entries currently allocated. */
  int capacity;
};

/**
 * @brief Grow a thread's ind scratch buffer if it is not already big enough.
 *
 * @param scratch The thread's #space_split_ind_scratch.
 * @param needed The number of entries required.
 */
static void space_split_ensure_ind_scratch(
    struct space_split_ind_scratch *scratch, const int needed) {

  if (needed <= scratch->capacity) return;

  if (scratch->ind != NULL) swift_free("ind_scratch", scratch->ind);
  if (swift_memalign("ind_scratch", (void **)&scratch->ind,
                     SWIFT_STRUCT_ALIGNMENT, sizeof(int) * needed) != 0)
    error("Failed to allocate ind scratch buffer.");
  scratch->capacity = needed;
}

/**
 * @brief Extra data for #space_split_mapper().
 *
 * The sorting buffers are allocated once, up front, for the complete local
 * particle arrays by #space_split() rather than per top-level cell (see
 * #space_split_fill_buffers()). This bundles the #space together with the
 * base pointer of each of those buffers and one #space_split_ind_scratch
 * and #space_split_particle_scratch per worker thread, so every worker can
 * find its own top-level cells' slices and its own reusable scratch space.
 */
struct space_split_mapper_data {

  /*! The #space being split. */
  struct space *s;

  /*! Base pointer of the space-wide hydro sorting buffer. */
  struct cell_buff *buff;

  /*! Base pointer of the space-wide gravity sorting buffer. */
  struct cell_buff *gbuff;

  /*! Base pointer of the space-wide star sorting buffer. */
  struct cell_buff *sbuff;

  /*! Base pointer of the space-wide black hole sorting buffer. */
  struct cell_buff *bbuff;

  /*! Base pointer of the space-wide sink sorting buffer. */
  struct cell_buff *sink_buff;

  /*! Per-thread #cell_split() ind scratch buffers, one per worker thread. */
  struct space_split_ind_scratch *ind_scratch;

  /*! Per-thread particle scratch buffers, one per worker thread. */
  struct space_split_particle_scratch *particle_scratch;
};

/**
 * @brief #threadpool mapper function to split cells if they contain
 *        too many particles.
 *
 * For each top-level cell, this first fills that cell's slice of the
 * space-wide sorting buffers (see #space_split_fill_buffers()), ensures
 * this thread's #cell_split() and particle scratch buffers are big enough
 * for it, then builds the cell hierarchy, sorts particles into it and
 * gathers each leaf's particles into the particle scratch buffers as it is
 * found (see #space_split_move_leaf()). Once that whole subtree has been
 * gathered, this copies each family's scratch data straight back into the
 * real arrays with a single bulk memcpy. It does not compute any per-cell
 * statistics; see #space_split_aggregate_mapper() for the pass that
 * finalises leaves and derives cell statistics and multipoles once the
 * whole hierarchy has been built.
 *
 * @param map_data Pointer towards the top-cells.
 * @param num_cells The number of cells to treat.
 * @param extra_data Pointer to a #space_split_mapper_data.
 */
static void space_split_mapper(void *map_data, int num_cells,
                               void *extra_data) {

  /* Unpack the inputs. */
  struct space_split_mapper_data *data =
      (struct space_split_mapper_data *)extra_data;
  struct space *s = data->s;
  struct cell *cells_top = s->cells_top;
  int *local_cells_with_particles = (int *)map_data;

  /* Threadpool id of current thread. */
  short int tpid = threadpool_gettid();
  struct space_split_particle_scratch *scratch = &data->particle_scratch[tpid];

  /* Loop over the non-empty cells */
  for (int ind = 0; ind < num_cells; ind++) {

    /* Get this cell. */
    struct cell *c = &cells_top[local_cells_with_particles[ind]];

    /* A top-level cell's particles occupy a contiguous slice of the
     * complete local particle arrays, at the same offset used elsewhere
     * (e.g. #cell_split()). Find this cell's slice of each space-wide
     * buffer and fill it from the cell's own particles. */
    const ptrdiff_t parts_offset = c->hydro.parts - s->parts;
    const ptrdiff_t gparts_offset = c->grav.parts - s->gparts;
    const ptrdiff_t sparts_offset = c->stars.parts - s->sparts;
    const ptrdiff_t bparts_offset = c->black_holes.parts - s->bparts;
    const ptrdiff_t sinks_offset = c->sinks.parts - s->sinks;
    struct cell_buff *buff = data->buff + parts_offset;
    struct cell_buff *gbuff = data->gbuff + gparts_offset;
    struct cell_buff *sbuff = data->sbuff + sparts_offset;
    struct cell_buff *bbuff = data->bbuff + bparts_offset;
    struct cell_buff *sink_buff = data->sink_buff + sinks_offset;
    space_split_fill_buffers(c, buff, gbuff, sbuff, bbuff, sink_buff,
                             parts_offset, gparts_offset, sparts_offset,
                             bparts_offset, sinks_offset);

    /* This thread's cell_split() scratch buffer only needs to be as big as
     * this top-level cell's own largest count -- every descendant's count
     * can only be smaller. */
    int max_count = c->hydro.count;
    max_count = max(max_count, c->grav.count);
    max_count = max(max_count, c->stars.count);
    max_count = max(max_count, c->black_holes.count);
    max_count = max(max_count, c->sinks.count);
    space_split_ensure_ind_scratch(&data->ind_scratch[tpid], max_count);

    /* This thread's particle scratch buffers need to hold this whole
     * top-level cell's particles, per family. */
    space_split_ensure_particle_scratch(scratch, c->hydro.count,
                                        c->grav.count, c->stars.count,
                                        c->black_holes.count,
                                        c->sinks.count);

    /* Split this cell recursively, gathering each leaf's particles into
     * scratch as it is found. */
    space_split_recursive(s, c, scratch->parts, scratch->xparts,
                          scratch->sparts, scratch->bparts, scratch->sinks,
                          scratch->gparts, data->ind_scratch[tpid].ind, buff,
                          sbuff, bbuff, gbuff, sink_buff, tpid);

    /* This cell's entire subtree has now been gathered into scratch, in
     * final order: copy each family back into the real arrays, skipping any
     * particle that was already in its final slot (source index == final
     * index). Such particles were never gathered into scratch (their slot
     * is garbage) and never need moving -- exactly the work master avoids
     * for particles that have not changed cell. */
    for (int k = 0; k < c->hydro.count; k++) {
      if (buff[k].part_ind == (size_t)(parts_offset + k)) continue;
      c->hydro.parts[k] = scratch->parts[k];
      c->hydro.xparts[k] = scratch->xparts[k];
    }
    for (int k = 0; k < c->stars.count; k++) {
      if (sbuff[k].part_ind == (size_t)(sparts_offset + k)) continue;
      c->stars.parts[k] = scratch->sparts[k];
    }
    for (int k = 0; k < c->black_holes.count; k++) {
      if (bbuff[k].part_ind == (size_t)(bparts_offset + k)) continue;
      c->black_holes.parts[k] = scratch->bparts[k];
    }
    for (int k = 0; k < c->sinks.count; k++) {
      if (sink_buff[k].part_ind == (size_t)(sinks_offset + k)) continue;
      c->sinks.parts[k] = scratch->sinks[k];
    }
    for (int k = 0; k < c->grav.count; k++) {
      if (gbuff[k].part_ind == (size_t)(gparts_offset + k)) continue;
      c->grav.parts[k] = scratch->gparts[k];
    }
  }
}

/**
 * @brief Extra data for #space_split_aggregate_mapper().
 */
struct space_split_aggregate_mapper_data {

  /*! The #space being aggregated. */
  struct space *s;

  /*! Per-thread particle scratch buffers, freed here (in parallel, on the
   * owning thread) now the split pass no longer needs them. */
  struct space_split_particle_scratch *particle_scratch;
};

/**
 * @brief #threadpool mapper function for finalising leaves and accumulating
 * cell properties after splitting.
 *
 * This function will loop over all top level cells and call
 * #space_split_aggregate_recursive on each one, which will recurse all the
 * way to the leaves, finalise them from their own particles, and then
 * accumulate cell properties and populate multipoles from the bottom up.
 *
 * Before that, it frees this worker thread's particle scratch buffers,
 * which the split pass is done with: doing it here runs the frees in
 * parallel across workers, each on the thread (and NUMA node) that
 * allocated and touched its own buffers, instead of serially on the main
 * thread. See #space_split_free_particle_scratch() for the idempotency
 * that makes this safe across multiple chunks per thread.
 *
 * @param map_data Pointer to the start of a chunk of cell indices.
 * @param num_cells Number of indices in this chunk.
 * @param extra_data Pointer to a #space_split_aggregate_mapper_data.
 */
static void space_split_aggregate_mapper(void *map_data, int num_cells,
                                         void *extra_data) {

  /* Unpack the inputs. */
  struct space_split_aggregate_mapper_data *data =
      (struct space_split_aggregate_mapper_data *)extra_data;
  struct space *s = data->s;
  struct cell *cells_top = s->cells_top;
  int *local_cells_with_particles = (int *)map_data;

  /* Free this thread's particle scratch, no longer needed after the split
   * pass -- in parallel, on the thread that owns it. */
  space_split_free_particle_scratch(
      &data->particle_scratch[threadpool_gettid()]);

  /* Initialise some global information about the top-level m-poles */
  float min_a_grav = FLT_MAX;
  float max_softening = 0.f;
  float max_mpole_power[SELF_GRAVITY_MULTIPOLE_ORDER + 1] = {0.f};

  /* Loop over the non-empty cells */
  for (int ind = 0; ind < num_cells; ind++) {

    /* Get this cell and aggregate its progeny's properties into it. */
    struct cell *c = &cells_top[local_cells_with_particles[ind]];
    space_split_aggregate_recursive(s, c);

    /* If we are running with self-gravity, collect the global min/max values
     * of the multipole properties. */
    if (s->with_self_gravity) {
      min_a_grav =
          min(min_a_grav, c->grav.multipole->m_pole.min_old_a_grav_norm);
      max_softening =
          max(max_softening, c->grav.multipole->m_pole.max_softening);
      for (int n = 0; n < SELF_GRAVITY_MULTIPOLE_ORDER + 1; ++n)
        max_mpole_power[n] =
            max(max_mpole_power[n], c->grav.multipole->m_pole.power[n]);
    }

    /* The deepest cell in this top tree determines the contribution to
     * the global max depth. */
    atomic_max(&s->maxdepth, c->maxdepth);
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* All cells and particles should have consistent h_max values. */
  for (int ind = 0; ind < num_cells; ind++) {
    int depth = 0;
    const struct cell *c = &cells_top[local_cells_with_particles[ind]];
    if (!checkCellhdxmax(c, &depth)) message("    at cell depth %d", depth);
  }
#endif

  /* Update the global min/max values of the multipole properties. */
  atomic_min_f(&s->min_a_grav, min_a_grav);
  atomic_max_f(&s->max_softening, max_softening);
  for (int n = 0; n < SELF_GRAVITY_MULTIPOLE_ORDER + 1; ++n)
    atomic_max_f(&s->max_mpole_power[n], max_mpole_power[n]);
}

/**
 * @brief Allocate the space-wide sorting buffers for the complete local
 * particle arrays.
 *
 * Called once by #space_split(), up front, before the split pass. Each
 * top-level cell later fills and uses its own slice of these buffers (see
 * #space_split_mapper() and #space_split_fill_buffers()); this only
 * allocates them, sized to the space's complete local particle arrays, for
 * each non-empty particle family.
 *
 * @param s The #space whose local particle arrays size the buffers.
 * @param buff Pointer to the hydro buffer pointer, allocated in place.
 * @param gbuff Pointer to the gravity buffer pointer, allocated in place.
 * @param sbuff Pointer to the star buffer pointer, allocated in place.
 * @param bbuff Pointer to the black hole buffer pointer, allocated in place.
 * @param sink_buff Pointer to the sink buffer pointer, allocated in place.
 */
static void space_split_allocate_buffers(struct space *s,
                                         struct cell_buff **buff,
                                         struct cell_buff **gbuff,
                                         struct cell_buff **sbuff,
                                         struct cell_buff **bbuff,
                                         struct cell_buff **sink_buff) {

  /* Only allocate a buffer for a family if the space holds any such
   * particles locally. */
  if (s->nr_parts > 0) {
    if (swift_memalign("tempbuff", (void **)buff, SWIFT_STRUCT_ALIGNMENT,
                       sizeof(struct cell_buff) * s->nr_parts) != 0)
      error("Failed to allocate temporary indices.");
  }
  if (s->nr_gparts > 0) {
    if (swift_memalign("tempgbuff", (void **)gbuff, SWIFT_STRUCT_ALIGNMENT,
                       sizeof(struct cell_buff) * s->nr_gparts) != 0)
      error("Failed to allocate temporary indices.");
  }
  if (s->nr_sparts > 0) {
    if (swift_memalign("tempsbuff", (void **)sbuff, SWIFT_STRUCT_ALIGNMENT,
                       sizeof(struct cell_buff) * s->nr_sparts) != 0)
      error("Failed to allocate temporary indices.");
  }
  if (s->nr_bparts > 0) {
    if (swift_memalign("tempbbuff", (void **)bbuff, SWIFT_STRUCT_ALIGNMENT,
                       sizeof(struct cell_buff) * s->nr_bparts) != 0)
      error("Failed to allocate temporary indices.");
  }
  if (s->nr_sinks > 0) {
    if (swift_memalign("temp_sink_buff", (void **)sink_buff,
                       SWIFT_STRUCT_ALIGNMENT,
                       sizeof(struct cell_buff) * s->nr_sinks) != 0)
      error("Failed to allocate temporary indices.");
  }
}

/**
 * @brief Split particles between cells of a hierarchy.
 *
 * This is done in parallel using threads in the #threadpool.
 * Only do this for the local non-empty top-level cells.
 *
 * The work is split into two passes, both run in parallel via the
 * #threadpool:
 *
 * 1. Split pass: #space_split_mapper() recursively builds the cell
 *    hierarchy and sorts each cell's slice of the space-wide sorting
 *    buffers into its progeny (see #cell_split()), gathering each leaf's
 *    particles into this thread's particle scratch buffers as it is found
 *    (see #space_split_move_leaf()), then, once a top-level cell's whole
 *    subtree has been gathered, copies each family's scratch data back
 *    into the real arrays with a single bulk memcpy. It does not compute
 *    any per-cell statistics.
 * 2. Aggregation pass: #space_split_aggregate_mapper() walks the
 *    now-complete hierarchy bottom-up, finalising leaves from their own
 *    (already moved) particles and deriving every other cell's statistics
 *    (h_max, time-step bounds, star formation history, maxdepth and, for
 *    self-gravity runs, the multipole) from its progeny.
 *
 * There is no array swap or pointer rebase between the two: every family
 * is copied back into the same array, at the same addresses, it always
 * occupied, so every cell pointer computed during the split pass is still
 * valid throughout the aggregation pass.
 *
 * @param s The #space.
 * @param verbose Are we talkative ?
 */
void space_split(struct space *s, int verbose) {

  const ticks tic = getticks();

  s->min_a_grav = FLT_MAX;
  s->max_softening = 0.f;
  bzero(s->max_mpole_power, (SELF_GRAVITY_MULTIPOLE_ORDER + 1) * sizeof(float));

  /* Allocate the sorting buffers for the complete local particle arrays;
   * each top-level cell will fill and use its own slice (see
   * #space_split_mapper()). */
  const ticks tic_alloc = getticks();
  struct cell_buff *buff = NULL, *gbuff = NULL, *sbuff = NULL, *bbuff = NULL,
                   *sink_buff = NULL;
  space_split_allocate_buffers(s, &buff, &gbuff, &sbuff, &bbuff, &sink_buff);

  /* One #cell_split() ind scratch buffer and one particle scratch buffer
   * per worker thread, both grown lazily as #space_split_mapper()
   * encounters top-level cells with larger counts -- the particle scratch
   * buffers only ever need to be as big as a single top-level cell's
   * counts, not the space's complete local particle arrays. */
  const int nr_threads = s->e->threadpool.num_threads;
  struct space_split_ind_scratch *ind_scratch =
      (struct space_split_ind_scratch *)swift_malloc(
          "ind_scratch_array",
          sizeof(struct space_split_ind_scratch) * nr_threads);
  if (ind_scratch == NULL) error("Failed to allocate ind scratch array.");
  bzero(ind_scratch, sizeof(struct space_split_ind_scratch) * nr_threads);
  struct space_split_particle_scratch *particle_scratch =
      (struct space_split_particle_scratch *)swift_malloc(
          "particle_scratch_array",
          sizeof(struct space_split_particle_scratch) * nr_threads);
  if (particle_scratch == NULL)
    error("Failed to allocate particle scratch array.");
  bzero(particle_scratch,
        sizeof(struct space_split_particle_scratch) * nr_threads);
  if (verbose)
    message("Allocation: %.3f %s.", clocks_from_ticks(getticks() - tic_alloc),
            clocks_getunit());

  /* Split pass: fill the buffers, build the cell hierarchy, sort the
   * particles into it, gather each leaf's particles into scratch as it is
   * found and copy each top-level cell's worth back into place once its
   * whole subtree has been gathered. */
  const ticks tic_split = getticks();
  struct space_split_mapper_data split_data = {
      s,     buff,        gbuff,          sbuff,
      bbuff, sink_buff,   ind_scratch,    particle_scratch};
  threadpool_map(&s->e->threadpool, space_split_mapper,
                 s->local_cells_with_particles_top,
                 s->nr_local_cells_with_particles, sizeof(int),
                 threadpool_auto_chunk_size, &split_data);
  if (verbose)
    message("Split pass: %.3f %s.", clocks_from_ticks(getticks() - tic_split),
            clocks_getunit());

#ifdef SWIFT_DEBUG_CHECKS
  part_verify_links(s->parts, s->gparts, s->sinks, s->sparts, s->bparts,
                    s->nr_parts, s->nr_gparts, s->nr_sinks, s->nr_sparts,
                    s->nr_bparts, /*verbose=*/0);
#endif

  /* The sorting buffers and the ind scratch are only needed for the split
   * pass above. The per-thread particle scratch buffers are freed in the
   * aggregation pass below, in parallel on their owning threads, so only
   * their (tiny) containing array is freed here -- after that pass. */
  const ticks tic_cleanup = getticks();
  if (buff != NULL) swift_free("tempbuff", buff);
  if (gbuff != NULL) swift_free("tempgbuff", gbuff);
  if (sbuff != NULL) swift_free("tempsbuff", sbuff);
  if (bbuff != NULL) swift_free("tempbbuff", bbuff);
  if (sink_buff != NULL) swift_free("temp_sink_buff", sink_buff);
  for (int t = 0; t < nr_threads; t++) {
    if (ind_scratch[t].ind != NULL) swift_free("ind_scratch", ind_scratch[t].ind);
  }
  swift_free("ind_scratch_array", ind_scratch);
  if (verbose)
    message("Buffer cleanup: %.3f %s.",
            clocks_from_ticks(getticks() - tic_cleanup), clocks_getunit());

  /* Aggregation pass: finalise leaves and derive cell statistics and
   * multipoles bottom-up (and free each thread's particle scratch in
   * parallel on the way). No rebase needed -- see #space_split()'s
   * documentation. */
  const ticks tic_aggregate = getticks();
  struct space_split_aggregate_mapper_data aggregate_data = {s,
                                                             particle_scratch};
  threadpool_map(&s->e->threadpool, space_split_aggregate_mapper,
                 s->local_cells_with_particles_top,
                 s->nr_local_cells_with_particles, sizeof(int),
                 threadpool_auto_chunk_size, &aggregate_data);

  /* Straggler sweep: free the scratch of any thread the aggregation pass
   * never handed a chunk to (a no-op for every thread that was), then the
   * containing array itself. */
  for (int t = 0; t < nr_threads; t++)
    space_split_free_particle_scratch(&particle_scratch[t]);
  swift_free("particle_scratch_array", particle_scratch);
  if (verbose)
    message("Aggregate pass: %.3f %s.",
            clocks_from_ticks(getticks() - tic_aggregate), clocks_getunit());

  if (verbose) {
    message("Max tree depth after split: %d", s->maxdepth);
    message("Have %d cells including subcells (cell footprint: %zd MB)",
            s->tot_cells, s->tot_cells * sizeof(struct cell) / (1024 * 1024));
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
  }
}
