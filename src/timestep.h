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
 * @param time_base_inv The inverse of the system's minimal time-step.
 */
__attribute__((always_inline)) INLINE static integertime_t
make_integer_timestep(float new_dt, timebin_t old_bin, integertime_t ti_current,
                      double time_base_inv) {

  /* Convert to integer time */
  integertime_t new_dti = (integertime_t)(new_dt * time_base_inv);

  /* Current time-step */
  integertime_t current_dti = get_integer_timestep(old_bin);
  integertime_t ti_end = get_integer_time_end(ti_current, old_bin);

  /* Limit timestep increase */
  if (old_bin > 0) new_dti = min(new_dti, 2 * current_dti);

  /* Put this timestep on the time line */
  integertime_t dti_timeline = max_nr_timesteps;
  while (new_dti < dti_timeline) dti_timeline /= ((integertime_t)2);
  new_dti = dti_timeline;

  /* Make sure we are allowed to increase the timestep size */
  if (new_dti > current_dti) {
    if ((max_nr_timesteps - ti_end) % new_dti > 0) new_dti = current_dti;
  }

#ifdef SWIFT_DEBUG_CHECKS
  if (new_dti == 0) error("Computed an integer time-step of size 0");
#endif

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

  float new_dt_self = FLT_MAX, new_dt_ext = FLT_MAX;

  if (e->policy & engine_policy_external_gravity)
    new_dt_ext = external_gravity_timestep(e->time, e->external_potential,
                                           e->physical_constants, gp);

  const float a_hydro[3] = {0.f, 0.f, 0.f};
  if (e->policy & engine_policy_self_gravity)
    new_dt_ext = gravity_compute_timestep_self(
        gp, a_hydro, e->gravity_properties, e->cosmology);

  /* Take the minimum of all */
  float new_dt = min(new_dt_self, new_dt_ext);

  /* Apply the maximal displacement constraint (FLT_MAX  if non-cosmological)*/
  new_dt = min(new_dt, e->dt_max_RMS_displacement);

  /* Apply cosmology correction (This is 1 if non-cosmological) */
  new_dt *= e->cosmology->time_step_factor;

  /* Limit timestep within the allowed range */
  new_dt = min(new_dt, e->dt_max);
  if (new_dt < e->dt_min)
    error("gpart (id=%lld) wants a time-step (%e) below dt_min (%e)",
          gp->id_or_neg_offset, new_dt, e->dt_min);

  /* Convert to integer time */
  const integertime_t new_dti = make_integer_timestep(
      new_dt, gp->time_bin, e->ti_current, e->time_base_inv);

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
  const float new_dt_hydro =
      hydro_compute_timestep(p, xp, e->hydro_properties, e->cosmology);

  /* Compute the next timestep (cooling condition) */
  float new_dt_cooling = FLT_MAX;
  if (e->policy & engine_policy_cooling)
    new_dt_cooling =
        cooling_timestep(e->cooling_func, e->physical_constants, e->cosmology,
                         e->internal_units, e->hydro_properties, p, xp);

  /* Compute the next timestep (gravity condition) */
  float new_dt_grav = FLT_MAX, new_dt_self_grav = FLT_MAX,
        new_dt_ext_grav = FLT_MAX;
  if (p->gpart != NULL) {

    if (e->policy & engine_policy_external_gravity)
      new_dt_ext_grav = external_gravity_timestep(
          e->time, e->external_potential, e->physical_constants, p->gpart);

    if (e->policy & engine_policy_self_gravity)
      new_dt_self_grav = gravity_compute_timestep_self(
          p->gpart, p->a_hydro, e->gravity_properties, e->cosmology);

    new_dt_grav = min(new_dt_self_grav, new_dt_ext_grav);
  }

  /* Final time-step is minimum of hydro and gravity */
  float new_dt = min3(new_dt_hydro, new_dt_cooling, new_dt_grav);

  /* Limit change in smoothing length */
  const float dt_h_change =
      (p->force.h_dt != 0.0f)
          ? fabsf(e->hydro_properties->log_max_h_change * p->h / p->force.h_dt)
          : FLT_MAX;

  new_dt = min(new_dt, dt_h_change);

  /* Apply the maximal displacement constraint (FLT_MAX if non-cosmological)*/
  new_dt = min(new_dt, e->dt_max_RMS_displacement);

  /* Apply cosmology correction (This is 1 if non-cosmological) */
  new_dt *= e->cosmology->time_step_factor;

  /* Limit timestep within the allowed range */
  new_dt = min(new_dt, e->dt_max);

  if (new_dt < e->dt_min)
    error("part (id=%lld) wants a time-step (%e) below dt_min (%e)", p->id,
          new_dt, e->dt_min);

  /* Convert to integer time */
  const integertime_t new_dti = make_integer_timestep(
      new_dt, p->time_bin, e->ti_current, e->time_base_inv);

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

  /* Stellar time-step */
  float new_dt_stars = stars_compute_timestep(sp);

  /* Gravity time-step */
  float new_dt_self = FLT_MAX, new_dt_ext = FLT_MAX;

  if (e->policy & engine_policy_external_gravity)
    new_dt_ext = external_gravity_timestep(e->time, e->external_potential,
                                           e->physical_constants, sp->gpart);

  const float a_hydro[3] = {0.f, 0.f, 0.f};
  if (e->policy & engine_policy_self_gravity)
    new_dt_self = gravity_compute_timestep_self(
        sp->gpart, a_hydro, e->gravity_properties, e->cosmology);

  /* Limit change in smoothing length */
  const float dt_h_change = (sp->feedback.h_dt != 0.0f)
                                ? fabsf(e->stars_properties->log_max_h_change *
                                        sp->h / sp->feedback.h_dt)
                                : FLT_MAX;

  /* Take the minimum of all */
  float new_dt = min4(new_dt_stars, new_dt_self, new_dt_ext, dt_h_change);

  /* Apply the maximal displacement constraint (FLT_MAX  if non-cosmological)*/
  new_dt = min(new_dt, e->dt_max_RMS_displacement);

  /* Apply cosmology correction (This is 1 if non-cosmological) */
  new_dt *= e->cosmology->time_step_factor;

  /* Limit timestep within the allowed range */
  new_dt = min(new_dt, e->dt_max);
  if (new_dt < e->dt_min) {
    error("spart (id=%lld) wants a time-step (%e) below dt_min (%e)", sp->id,
          new_dt, e->dt_min);
  }

  /* Convert to integer time */
  const integertime_t new_dti = make_integer_timestep(
      new_dt, sp->time_bin, e->ti_current, e->time_base_inv);

  return new_dti;
}

#endif /* SWIFT_TIMESTEP_H */
