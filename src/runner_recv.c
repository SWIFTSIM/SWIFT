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

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "runner.h"

/* Local headers. */
#include "engine.h"
#include "timers.h"

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
  const timebin_t max_active_bin = r->e->max_active_bin;

  TIMER_TIC;

  integertime_t ti_hydro_end_min = max_nr_timesteps;
  integertime_t ti_hydro_end_max = 0;
  timebin_t time_bin_min = num_time_bins;
  timebin_t time_bin_max = 0;
  float h_max = 0.f;
  float h_max_active = 0.f;

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
      if (parts[k].time_bin <= max_active_bin)
        h_max_active = max(h_max_active, parts[k].h);
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
        h_max_active = max(h_max_active, c->progeny[k]->hydro.h_max_active);
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (!(r->e->policy & engine_policy_timestep_sync) &&
      !(r->e->policy & engine_policy_timestep_limiter) &&
      ti_hydro_end_min < ti_current)
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
  c->hydro.h_max_active = h_max_active;

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
  if (!(r->e->policy & engine_policy_timestep_sync) &&
      !(r->e->policy & engine_policy_timestep_limiter) &&
      ti_gravity_end_min < ti_current)
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
  const timebin_t max_active_bin = r->e->max_active_bin;

  TIMER_TIC;

  integertime_t ti_stars_end_min = max_nr_timesteps;
  integertime_t ti_stars_end_max = 0;
  timebin_t time_bin_min = num_time_bins;
  timebin_t time_bin_max = 0;
  float h_max = 0.f;
  float h_max_active = 0.f;

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
      sparts[k].num_ngb_feedback = 0;
#endif
      if (sparts[k].time_bin == time_bin_inhibited) continue;
      time_bin_min = min(time_bin_min, sparts[k].time_bin);
      time_bin_max = max(time_bin_max, sparts[k].time_bin);
      h_max = max(h_max, sparts[k].h);
      if (sparts[k].time_bin <= max_active_bin)
        h_max_active = max(h_max_active, sparts[k].h);
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
            max(ti_stars_end_max, c->progeny[k]->stars.ti_end_max);
        h_max = max(h_max, c->progeny[k]->stars.h_max);
        h_max_active = max(h_max_active, c->progeny[k]->stars.h_max_active);
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_stars_end_min < ti_current &&
      !(r->e->policy & engine_policy_star_formation))
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
  c->stars.h_max_active = h_max_active;

  if (timer) TIMER_TOC(timer_dorecv_spart);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}

/**
 * @brief Construct the cell properties from the received #bpart.
 *
 * Note that we do not need to clear the sorts since we do not sort
 * the black holes.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param clear_sorts Should we clear the sort flag and hence trigger a sort ?
 * @param timer Are we timing this ?
 */
void runner_do_recv_bpart(struct runner *r, struct cell *c, int clear_sorts,
                          int timer) {

#ifdef WITH_MPI

  struct bpart *restrict bparts = c->black_holes.parts;
  const size_t nr_bparts = c->black_holes.count;
  const integertime_t ti_current = r->e->ti_current;
  const timebin_t max_active_bin = r->e->max_active_bin;

  TIMER_TIC;

  integertime_t ti_black_holes_end_min = max_nr_timesteps;
  integertime_t ti_black_holes_end_max = 0;
  timebin_t time_bin_min = num_time_bins;
  timebin_t time_bin_max = 0;
  float h_max = 0.f;
  float h_max_active = 0.f;

#ifdef SWIFT_DEBUG_CHECKS
  if (c->nodeID == engine_rank) error("Updating a local cell!");
#endif

  /* If this cell is a leaf, collect the particle data. */
  if (!c->split) {

    /* Collect everything... */
    for (size_t k = 0; k < nr_bparts; k++) {
#ifdef DEBUG_INTERACTIONS_BLACK_HOLES
      bparts[k].num_ngb_force = 0;
#endif

      if (bparts[k].time_bin == time_bin_inhibited) continue;
      time_bin_min = min(time_bin_min, bparts[k].time_bin);
      time_bin_max = max(time_bin_max, bparts[k].time_bin);
      h_max = max(h_max, bparts[k].h);
      if (bparts[k].time_bin <= max_active_bin)
        h_max_active = max(h_max_active, bparts[k].h);
    }

    /* Convert into a time */
    ti_black_holes_end_min = get_integer_time_end(ti_current, time_bin_min);
    ti_black_holes_end_max = get_integer_time_end(ti_current, time_bin_max);
  }

  /* Otherwise, recurse and collect. */
  else {
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL && c->progeny[k]->black_holes.count > 0) {
        runner_do_recv_bpart(r, c->progeny[k], clear_sorts, 0);
        ti_black_holes_end_min =
            min(ti_black_holes_end_min, c->progeny[k]->black_holes.ti_end_min);
        ti_black_holes_end_max =
            max(ti_black_holes_end_max, c->progeny[k]->black_holes.ti_end_max);
        h_max = max(h_max, c->progeny[k]->black_holes.h_max);
        h_max_active =
            max(h_max_active, c->progeny[k]->black_holes.h_max_active);
      }
    }
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (ti_black_holes_end_min < ti_current)
    error(
        "Received a cell at an incorrect time c->ti_end_min=%lld, "
        "e->ti_current=%lld.",
        ti_black_holes_end_min, ti_current);
#endif

  /* ... and store. */
  // c->grav.ti_end_min = ti_gravity_end_min;
  // c->grav.ti_end_max = ti_gravity_end_max;
  c->black_holes.ti_old_part = ti_current;
  c->black_holes.h_max = h_max;
  c->black_holes.h_max_active = h_max_active;

  if (timer) TIMER_TOC(timer_dorecv_bpart);

#else
  error("SWIFT was not compiled with MPI support.");
#endif
}
