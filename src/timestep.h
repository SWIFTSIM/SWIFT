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
#include "const.h"
#include "cooling.h"
#include "debug.h"
/**
 * @brief Compute a valid integer time-step form a given time-step
 *
 * @param new_dt The time-step to convert.
 * @param ti_begin The (integer) start of the previous time-step.
 * @param ti_end The (integer) end of the previous time-step.
 * @param timeBase_inv The inverse of the system's minimal time-step.
 */
__attribute__((always_inline)) INLINE static int get_integer_timestep(
    float new_dt, int ti_begin, int ti_end, double timeBase_inv) {

  /* Convert to integer time */
  int new_dti = (int)(new_dt * timeBase_inv);

  /* Recover the current timestep */
  const int current_dti = ti_end - ti_begin;

  /* Limit timestep increase */
  if (current_dti > 0) new_dti = min(new_dti, 2 * current_dti);

  /* Put this timestep on the time line */
  int dti_timeline = max_nr_timesteps;
  while (new_dti < dti_timeline) dti_timeline /= 2;
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
__attribute__((always_inline)) INLINE static int get_gpart_timestep(
    const struct gpart *restrict gp, const struct engine *restrict e) {

  const float new_dt_external = external_gravity_timestep(
      e->time, e->external_potential, e->physical_constants, gp);

  /* const float new_dt_self = */
  /*     gravity_compute_timestep_self(e->physical_constants, gp); */
  const float new_dt_self = FLT_MAX;  // MATTHIEU

  float new_dt = min(new_dt_external, new_dt_self);

  /* Limit timestep within the allowed range */
  new_dt = min(new_dt, e->dt_max);
  new_dt = max(new_dt, e->dt_min);

  /* Convert to integer time */
  const int new_dti =
      get_integer_timestep(new_dt, gp->ti_begin, gp->ti_end, e->timeBase_inv);

  return new_dti;
}

/**
 * @brief Compute the new (integer) time-step of a given #part
 *
 * @param p The #part.
 * @param xp The #xpart partner of p.
 * @param e The #engine (used to get some constants).
 */
__attribute__((always_inline)) INLINE static int get_part_timestep(
    const struct part *restrict p, const struct xpart *restrict xp,
    const struct engine *restrict e) {

  /* Compute the next timestep (hydro condition) */
  const float new_dt_hydro = hydro_compute_timestep(p, xp, e->hydro_properties);

  /* Compute the next timestep (cooling condition) */
  float new_dt_cooling = FLT_MAX;
  if (e->policy & engine_policy_cooling)
    new_dt_cooling = cooling_timestep(e->cooling_func, e->physical_constants,
                                      e->internalUnits, p);

  /* Compute the next timestep (gravity condition) */
  float new_dt_grav = FLT_MAX;
  if (p->gpart != NULL) {

    const float new_dt_external = external_gravity_timestep(
        e->time, e->external_potential, e->physical_constants, p->gpart);

    /* const float new_dt_self = */
    /*     gravity_compute_timestep_self(e->physical_constants, p->gpart); */
    const float new_dt_self = FLT_MAX;  // MATTHIEU

    new_dt_grav = min(new_dt_external, new_dt_self);
  }

  /* Final time-step is minimum of hydro and gravity */
  float new_dt = min(min(new_dt_hydro, new_dt_cooling), new_dt_grav);

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
  const int new_dti =
      get_integer_timestep(new_dt, p->ti_begin, p->ti_end, e->timeBase_inv);

  return new_dti;
}

#endif /* SWIFT_TIMESTEP_H */
