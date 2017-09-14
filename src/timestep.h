/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_TIMESTEP_H
#define SWIFT_TIMESTEP_H

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "cooling.h"
#include "debug.h"
#include "timeline.h"

/**
 * @brief Compute a valid integer time-step form a given time-step
 *
 * @param new_dt The time-step to convert.
 * @param old_bin The old time bin.
 * @param ti_current The current time on the integer time-line.
 * @param timeBase_inv The inverse of the system's minimal time-step.
 */
__attribute__((always_inline)) INLINE static integertime_t
make_integer_timestep(float new_dt, timebin_t old_bin, integertime_t ti_current,
                      double timeBase_inv) {

  /* Convert to integer time */
  integertime_t new_dti = (integertime_t)(new_dt * timeBase_inv);

  /* Current time-step */
  integertime_t current_dti = get_integer_timestep(old_bin);
  integertime_t ti_end = get_integer_time_end(ti_current, old_bin);

  /* Limit timestep increase */
  if (old_bin > 0) new_dti = min(new_dti, 2 * current_dti);

  /* Put this timestep on the time line */
  integertime_t dti_timeline = max_nr_timesteps;
  while (new_dti < dti_timeline) dti_timeline /= 2LL;
  new_dti = dti_timeline;

  /* Make sure we are allowed to increase the timestep size */
  if (new_dti > current_dti) {
    if ((max_nr_timesteps - ti_end) % new_dti > 0) new_dti = current_dti;
  }
  return new_dti;
}

/**
 * @brief Compute the new (integer) time-step of a given #gpart
 *
 * @param gp The #gpart.
 * @param e The #engine (used to get some constants).
 */
__attribute__((always_inline)) INLINE static integertime_t get_gpart_timestep(
    const struct gpart *restrict gp, const struct engine *restrict e) {

  float new_dt = FLT_MAX;

  if (e->policy & engine_policy_external_gravity)
    new_dt =
        min(new_dt, external_gravity_timestep(e->time, e->external_potential,
                                              e->physical_constants, gp));

  if (e->policy & engine_policy_self_gravity)
    new_dt =
        min(new_dt, gravity_compute_timestep_self(gp, e->gravity_properties));

  /* Limit timestep within the allowed range */
  new_dt = min(new_dt, e->dt_max);
  new_dt = max(new_dt, e->dt_min);

  /* Convert to integer time */
  const integertime_t new_dti = make_integer_timestep(
      new_dt, gp->time_bin, e->ti_current, e->timeBase_inv);

  return new_dti;
}

/**
 * @brief Compute the new (integer) time-step of a given #part
 *
 * @param p The #part.
 * @param xp The #xpart partner of p.
 * @param e The #engine (used to get some constants).
 */
__attribute__((always_inline)) INLINE static integertime_t get_part_timestep(
    const struct part *restrict p, const struct xpart *restrict xp,
    const struct engine *restrict e) {

  /* Compute the next timestep (hydro condition) */
  const float new_dt_hydro = hydro_compute_timestep(p, xp, e->hydro_properties);

  /* Compute the next timestep (cooling condition) */
  float new_dt_cooling = FLT_MAX;
  if (e->policy & engine_policy_cooling)
    new_dt_cooling = cooling_timestep(e->cooling_func, e->physical_constants,
                                      e->internal_units, p);

  /* Compute the next timestep (gravity condition) */
  float new_dt_grav = FLT_MAX;
  if (p->gpart != NULL) {

    if (e->policy & engine_policy_external_gravity)
      new_dt_grav = min(new_dt_grav, external_gravity_timestep(
                                         e->time, e->external_potential,
                                         e->physical_constants, p->gpart));

    if (e->policy & engine_policy_self_gravity)
      new_dt_grav = min(new_dt_grav, gravity_compute_timestep_self(
                                         p->gpart, e->gravity_properties));
  }

  /* Final time-step is minimum of hydro and gravity */
  float new_dt = min3(new_dt_hydro, new_dt_cooling, new_dt_grav);

  /* Limit change in h */
  const float dt_h_change =
      (p->force.h_dt != 0.0f)
          ? fabsf(e->hydro_properties->log_max_h_change * p->h / p->force.h_dt)
          : FLT_MAX;

  new_dt = min(new_dt, dt_h_change);

  /* Limit timestep within the allowed range */
  new_dt = min(new_dt, e->dt_max);
  new_dt = max(new_dt, e->dt_min);

  /* Convert to integer time */
  const integertime_t new_dti = make_integer_timestep(
      new_dt, p->time_bin, e->ti_current, e->timeBase_inv);

  return new_dti;
}

/**
 * @brief Compute the new (integer) time-step of a given #spart
 *
 * @param sp The #spart.
 * @param e The #engine (used to get some constants).
 */
__attribute__((always_inline)) INLINE static integertime_t get_spart_timestep(
    const struct spart *restrict sp, const struct engine *restrict e) {

  float new_dt = star_compute_timestep(sp);

  if (e->policy & engine_policy_external_gravity)
    new_dt = min(new_dt,
                 external_gravity_timestep(e->time, e->external_potential,
                                           e->physical_constants, sp->gpart));

  if (e->policy & engine_policy_self_gravity)
    new_dt = min(new_dt, gravity_compute_timestep_self(sp->gpart,
                                                       e->gravity_properties));

  /* Limit timestep within the allowed range */
  new_dt = min(new_dt, e->dt_max);
  new_dt = max(new_dt, e->dt_min);

  /* Convert to integer time */
  const integertime_t new_dti = make_integer_timestep(
      new_dt, sp->time_bin, e->ti_current, e->timeBase_inv);

  return new_dti;
}

#endif /* SWIFT_TIMESTEP_H */
