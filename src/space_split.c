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
 * @brief Recursively split a cell.
 *
 * @param s The #space in which the cell lives.
 * @param c The #cell to split recursively.
 * @param buff A buffer for particle sorting, should be of size at least
 *        c->hydro.count or @c NULL.
 * @param sbuff A buffer for particle sorting, should be of size at least
 *        c->stars.count or @c NULL.
 * @param bbuff A buffer for particle sorting, should be of size at least
 *        c->black_holes.count or @c NULL.
 * @param gbuff A buffer for particle sorting, should be of size at least
 *        c->grav.count or @c NULL.
 * @param sink_buff A buffer for particle sorting, should be of size at least
 *        c->sinks.count or @c NULL.
 */
void space_split_recursive(struct space *s, struct cell *c,
                           struct cell_buff *restrict buff,
                           struct cell_buff *restrict sbuff,
                           struct cell_buff *restrict bbuff,
                           struct cell_buff *restrict gbuff,
                           struct cell_buff *restrict sink_buff,
                           const short int tpid) {

  const int count = c->hydro.count;
  const int gcount = c->grav.count;
  const int scount = c->stars.count;
  const int bcount = c->black_holes.count;
  const int sink_count = c->sinks.count;
  const int with_self_gravity = s->with_self_gravity;
  const int depth = c->depth;
  int maxdepth = 0;
  float h_max = 0.0f;
  float h_max_active = 0.0f;
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
  struct part *parts = c->hydro.parts;
  struct gpart *gparts = c->grav.parts;
  struct spart *sparts = c->stars.parts;
  struct bpart *bparts = c->black_holes.parts;
  struct xpart *xparts = c->hydro.xparts;
  struct sink *sinks = c->sinks.parts;
  struct engine *e = s->e;
  const integertime_t ti_current = e->ti_current;
  const int with_rt = e->policy & engine_policy_rt;

  /* Set the top level cell tpid. Doing it here ensures top level cells
   * have the same tpid as their progeny. */
  if (depth == 0) c->tpid = tpid;

  /* If the buff is NULL, allocate it, and remember to free it. */
  const int allocate_buffer = (buff == NULL && gbuff == NULL && sbuff == NULL &&
                               bbuff == NULL && sink_buff == NULL);
  if (allocate_buffer) {
    if (count > 0) {
      if (swift_memalign("tempbuff", (void **)&buff, SWIFT_STRUCT_ALIGNMENT,
                         sizeof(struct cell_buff) * count) != 0)
        error("Failed to allocate temporary indices.");
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
      }
    }
    if (gcount > 0) {
      if (swift_memalign("tempgbuff", (void **)&gbuff, SWIFT_STRUCT_ALIGNMENT,
                         sizeof(struct cell_buff) * gcount) != 0)
        error("Failed to allocate temporary indices.");
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
      }
    }
    if (scount > 0) {
      if (swift_memalign("tempsbuff", (void **)&sbuff, SWIFT_STRUCT_ALIGNMENT,
                         sizeof(struct cell_buff) * scount) != 0)
        error("Failed to allocate temporary indices.");
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
      }
    }
    if (bcount > 0) {
      if (swift_memalign("tempbbuff", (void **)&bbuff, SWIFT_STRUCT_ALIGNMENT,
                         sizeof(struct cell_buff) * bcount) != 0)
        error("Failed to allocate temporary indices.");
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
      }
    }
    if (sink_count > 0) {
      if (swift_memalign("temp_sink_buff", (void **)&sink_buff,
                         SWIFT_STRUCT_ALIGNMENT,
                         sizeof(struct cell_buff) * sink_count) != 0)
        error("Failed to allocate temporary indices.");
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
      }
    }
  }

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
      cp->sinks.r_cut_max = 0.f;
      cp->sinks.r_cut_max_active = 0.f;
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
    cell_split(c, c->hydro.parts - s->parts, c->stars.parts - s->sparts,
               c->black_holes.parts - s->bparts, c->sinks.parts - s->sinks,
               buff, sbuff, bbuff, gbuff, sink_buff);

    /* Buffers for the progenitors */
    struct cell_buff *progeny_buff = buff, *progeny_gbuff = gbuff,
                     *progeny_sbuff = sbuff, *progeny_bbuff = bbuff,
                     *progeny_sink_buff = sink_buff;

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
        space_split_recursive(s, cp, progeny_buff, progeny_sbuff, progeny_bbuff,
                              progeny_gbuff, progeny_sink_buff, tpid);

        /* Update the pointers in the buffers */
        progeny_buff += cp->hydro.count;
        progeny_gbuff += cp->grav.count;
        progeny_sbuff += cp->stars.count;
        progeny_bbuff += cp->black_holes.count;
        progeny_sink_buff += cp->sinks.count;

        /* Update the cell-wide properties */
        h_max = max(h_max, cp->hydro.h_max);
        h_max_active = max(h_max_active, cp->hydro.h_max_active);
        stars_h_max = max(stars_h_max, cp->stars.h_max);
        stars_h_max_active = max(stars_h_max_active, cp->stars.h_max_active);
        black_holes_h_max = max(black_holes_h_max, cp->black_holes.h_max);
        black_holes_h_max_active =
            max(black_holes_h_max_active, cp->black_holes.h_max_active);
        sinks_h_max = max(sinks_h_max, cp->sinks.r_cut_max);
        sinks_h_max_active =
            max(sinks_h_max_active, cp->sinks.r_cut_max_active);

        ti_hydro_end_min = min(ti_hydro_end_min, cp->hydro.ti_end_min);
        ti_hydro_beg_max = max(ti_hydro_beg_max, cp->hydro.ti_beg_max);
        ti_rt_end_min = min(ti_rt_end_min, cp->rt.ti_rt_end_min);
        ti_rt_beg_max = max(ti_rt_beg_max, cp->rt.ti_rt_beg_max);
        ti_rt_min_step_size =
            min(ti_rt_min_step_size, cp->rt.ti_rt_min_step_size);
        ti_gravity_end_min = min(ti_gravity_end_min, cp->grav.ti_end_min);
        ti_gravity_beg_max = max(ti_gravity_beg_max, cp->grav.ti_beg_max);
        ti_stars_end_min = min(ti_stars_end_min, cp->stars.ti_end_min);
        ti_stars_beg_max = max(ti_stars_beg_max, cp->stars.ti_beg_max);
        ti_sinks_end_min = min(ti_sinks_end_min, cp->sinks.ti_end_min);
        ti_sinks_beg_max = max(ti_sinks_beg_max, cp->sinks.ti_beg_max);
        ti_black_holes_end_min =
            min(ti_black_holes_end_min, cp->black_holes.ti_end_min);
        ti_black_holes_beg_max =
            max(ti_black_holes_beg_max, cp->black_holes.ti_beg_max);

        star_formation_logger_add(&c->stars.sfh, &cp->stars.sfh);

        /* Increase the depth */
        maxdepth = max(maxdepth, cp->maxdepth);
      }
    }

    /* Deal with the multipole */
    if (s->with_self_gravity) {

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
          gravity_M2M(&temp, m, c->grav.multipole->CoM,
                      cp->grav.multipole->CoM);
          gravity_multipole_add(&c->grav.multipole->m_pole, &temp);

          /* Upper limit of max CoM<->gpart distance */
          const double dx =
              c->grav.multipole->CoM[0] - cp->grav.multipole->CoM[0];
          const double dy =
              c->grav.multipole->CoM[1] - cp->grav.multipole->CoM[1];
          const double dz =
              c->grav.multipole->CoM[2] - cp->grav.multipole->CoM[2];
          const double r2 = dx * dx + dy * dy + dz * dz;
          r_max = max(r_max, cp->grav.multipole->r_max + sqrt(r2));
        }
      }

      /* Alternative upper limit of max CoM<->gpart distance */
      const double dx =
          c->grav.multipole->CoM[0] > c->loc[0] + c->width[0] / 2.
              ? c->grav.multipole->CoM[0] - c->loc[0]
              : c->loc[0] + c->width[0] - c->grav.multipole->CoM[0];
      const double dy =
          c->grav.multipole->CoM[1] > c->loc[1] + c->width[1] / 2.
              ? c->grav.multipole->CoM[1] - c->loc[1]
              : c->loc[1] + c->width[1] - c->grav.multipole->CoM[1];
      const double dz =
          c->grav.multipole->CoM[2] > c->loc[2] + c->width[2] / 2.
              ? c->grav.multipole->CoM[2] - c->loc[2]
              : c->loc[2] + c->width[2] - c->grav.multipole->CoM[2];

      /* Take minimum of both limits */
      c->grav.multipole->r_max = min(r_max, sqrt(dx * dx + dy * dy + dz * dz));

      /* Store the value at rebuild time */
      c->grav.multipole->r_max_rebuild = c->grav.multipole->r_max;
      c->grav.multipole->CoM_rebuild[0] = c->grav.multipole->CoM[0];
      c->grav.multipole->CoM_rebuild[1] = c->grav.multipole->CoM[1];
      c->grav.multipole->CoM_rebuild[2] = c->grav.multipole->CoM[2];

      /* Compute the multipole power */
      gravity_multipole_compute_power(&c->grav.multipole->m_pole);

    } /* Deal with gravity */
  }   /* Split or let it be? */

  /* Otherwise, collect the data from the particles this cell. */
  else {

    /* Clear the progeny. */
    bzero(c->progeny, sizeof(struct cell *) * 8);
    c->split = 0;
    maxdepth = c->depth;

    ti_hydro_end_min = max_nr_timesteps;
    ti_hydro_beg_max = 0;

    ti_gravity_end_min = max_nr_timesteps;
    ti_gravity_beg_max = 0;

    ti_stars_end_min = max_nr_timesteps;
    ti_stars_beg_max = 0;

    ti_black_holes_end_min = max_nr_timesteps;
    ti_black_holes_beg_max = 0;

    ti_rt_end_min = max_nr_timesteps;
    ti_rt_beg_max = 0;
    ti_rt_min_step_size = max_nr_timesteps;

    /* parts: Get dt_min/dt_max and h_max. */
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

      h_max = max(h_max, parts[k].h);

      if (part_is_active(&parts[k], e))
        h_max_active = max(h_max_active, parts[k].h);

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

      stars_h_max = max(stars_h_max, sparts[k].h);

      if (spart_is_active(&sparts[k], e))
        stars_h_max_active = max(stars_h_max_active, sparts[k].h);

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

      sinks_h_max = max(sinks_h_max, sinks[k].r_cut);

      if (sink_is_active(&sinks[k], e))
        sinks_h_max_active = max(sinks_h_max_active, sinks[k].r_cut);

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

      black_holes_h_max = max(black_holes_h_max, bparts[k].h);

      if (bpart_is_active(&bparts[k], e))
        black_holes_h_max_active = max(black_holes_h_max_active, bparts[k].h);

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
    }
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
  c->sinks.r_cut_max = sinks_h_max;
  c->sinks.r_cut_max_active = sinks_h_max_active;
  c->black_holes.ti_end_min = ti_black_holes_end_min;
  c->black_holes.ti_beg_max = ti_black_holes_beg_max;
  c->black_holes.h_max = black_holes_h_max;
  c->black_holes.h_max_active = black_holes_h_max_active;
  c->maxdepth = maxdepth;

  /* No runner owns this cell yet. We assign those during scheduling. */
  c->owner = -1;

  /* Store the global max depth */
  if (c->depth == 0) atomic_max(&s->maxdepth, maxdepth);

  /* Clean up. */
  if (allocate_buffer) {
    if (buff != NULL) swift_free("tempbuff", buff);
    if (gbuff != NULL) swift_free("tempgbuff", gbuff);
    if (sbuff != NULL) swift_free("tempsbuff", sbuff);
    if (bbuff != NULL) swift_free("tempbbuff", bbuff);
    if (sink_buff != NULL) swift_free("temp_sink_buff", sink_buff);
  }
}

/**
 * @brief #threadpool mapper function to split cells if they contain
 *        too many particles.
 *
 * @param map_data Pointer towards the top-cells.
 * @param num_cells The number of cells to treat.
 * @param extra_data Pointers to the #space.
 */
void space_split_mapper(void *map_data, int num_cells, void *extra_data) {

  /* Unpack the inputs. */
  struct space *s = (struct space *)extra_data;
  struct cell *cells_top = s->cells_top;
  int *local_cells_with_particles = (int *)map_data;

  /* Collect some global information about the top-level m-poles */
  float min_a_grav = FLT_MAX;
  float max_softening = 0.f;
  float max_mpole_power[SELF_GRAVITY_MULTIPOLE_ORDER + 1] = {0.f};

  /* Threadpool id of current thread. */
  short int tpid = threadpool_gettid();

  /* Loop over the non-empty cells */
  for (int ind = 0; ind < num_cells; ind++) {
    struct cell *c = &cells_top[local_cells_with_particles[ind]];
    space_split_recursive(s, c, NULL, NULL, NULL, NULL, NULL, tpid);

    if (s->with_self_gravity) {
      min_a_grav =
          min(min_a_grav, c->grav.multipole->m_pole.min_old_a_grav_norm);
      max_softening =
          max(max_softening, c->grav.multipole->m_pole.max_softening);

      for (int n = 0; n < SELF_GRAVITY_MULTIPOLE_ORDER + 1; ++n)
        max_mpole_power[n] =
            max(max_mpole_power[n], c->grav.multipole->m_pole.power[n]);
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  /* All cells and particles should have consistent h_max values. */
  for (int ind = 0; ind < num_cells; ind++) {
    int depth = 0;
    const struct cell *c = &cells_top[local_cells_with_particles[ind]];
    if (!checkCellhdxmax(c, &depth)) message("    at cell depth %d", depth);
  }
#endif

  atomic_min_f(&s->min_a_grav, min_a_grav);
  atomic_max_f(&s->max_softening, max_softening);
  for (int n = 0; n < SELF_GRAVITY_MULTIPOLE_ORDER + 1; ++n)
    atomic_max_f(&s->max_mpole_power[n], max_mpole_power[n]);
}

/**
 * @brief Split particles between cells of a hierarchy.
 *
 * This is done in parallel using threads in the #threadpool.
 * Only do this for the local non-empty top-level cells.
 *
 * @param s The #space.
 * @param verbose Are we talkative ?
 */
void space_split(struct space *s, int verbose) {

  const ticks tic = getticks();

  s->min_a_grav = FLT_MAX;
  s->max_softening = 0.f;
  bzero(s->max_mpole_power, (SELF_GRAVITY_MULTIPOLE_ORDER + 1) * sizeof(float));

  threadpool_map(&s->e->threadpool, space_split_mapper,
                 s->local_cells_with_particles_top,
                 s->nr_local_cells_with_particles, sizeof(int),
                 threadpool_auto_chunk_size, s);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
