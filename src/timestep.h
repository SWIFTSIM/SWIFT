/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2016 Matthieu Schaller (schaller@strw.leidenuniv.nl)
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
#include <config.h>

/* Local headers. */
#include "cooling.h"
#include "debug.h"
#include "forcing.h"
#include "potential.h"
#include "rt.h"
#include "timeline.h"

/**
 * @brief Compute a valid integer time-step form a given time-step
 *
 * We consider the minimal time-bin of any neighbours and prevent particles
 * to differ from it by a fixed constant `time_bin_neighbour_max_delta_bin`.
 *
 * If min_ngb_bin is set to `num_time_bins`, then no limit from the neighbours
 * is imposed.
 *
 * @param new_dt The time-step to convert.
 * @param old_bin The old time bin.
 * @param min_ngb_bin Minimal time-bin of any neighbour of this particle.
 * @param ti_current The current time on the integer time-line.
 * @param time_base_inv The inverse of the system's minimal time-step.
 */
__attribute__((always_inline, const)) INLINE static integertime_t
make_integer_timestep(const float new_dt, const timebin_t old_bin,
                      const timebin_t min_ngb_bin,
                      const integertime_t ti_current,
                      const double time_base_inv) {

  /* Convert to integer time */
  integertime_t new_dti = (integertime_t)(new_dt * time_base_inv);

  /* Are we allowed to use this bin given the neighbours? */
  timebin_t new_bin = get_time_bin(new_dti);
  new_bin = min(new_bin, min_ngb_bin + time_bin_neighbour_max_delta_bin);
  new_dti = get_integer_timestep(new_bin);

  /* Current time-step */
  const integertime_t current_dti = get_integer_timestep(old_bin);
  const integertime_t ti_end = get_integer_time_end(ti_current, old_bin);

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

#ifdef SWIFT_DEBUG_CHECKS
  if (gp->time_bin == time_bin_not_created) {
    error("Trying to compute time step for an extra particle.");
  }
#endif

  float new_dt_self = FLT_MAX, new_dt_ext = FLT_MAX;

  if (e->policy & engine_policy_external_gravity)
    new_dt_ext = external_gravity_timestep(e->time, e->external_potential,
                                           e->physical_constants, gp);

  const float a_hydro[3] = {0.f, 0.f, 0.f};
  if (e->policy & engine_policy_self_gravity)
    new_dt_self = gravity_compute_timestep_self(
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
      new_dt, gp->time_bin, num_time_bins, e->ti_current, e->time_base_inv);

  return new_dti;
}

/**
 * @brief Compute the new (integer) time-step of a given #part
 *
 * @param p The #part.
 * @param xp The #xpart partner of p.
 * @param e The #engine (used to get some constants).
 * @param new_dti_rt The new radiation integer time step.
 */
__attribute__((always_inline)) INLINE static integertime_t get_part_timestep(
    const struct part *restrict p, const struct xpart *restrict xp,
    const struct engine *restrict e, const integertime_t new_dti_rt) {

  /* Compute the next timestep (hydro condition) */
  const float new_dt_hydro =
      hydro_compute_timestep(p, xp, e->hydro_properties, e->cosmology);

  /* Compute the next timestep (MHD condition) */
  const float new_dt_mhd =
      mhd_compute_timestep(p, xp, e->hydro_properties, e->cosmology);

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

  /* Compute the next timestep (forcing terms condition) */
  const float new_dt_forcing = forcing_terms_timestep(
      e->time, e->forcing_terms, e->physical_constants, p, xp);

  /* Compute the next timestep (chemistry condition, e.g. diffusion) */
  const float new_dt_chemistry =
      chemistry_timestep(e->physical_constants, e->cosmology, e->internal_units,
                         e->hydro_properties, e->chemistry, p);

  /* Take the minimum of all */
  float new_dt = min3(new_dt_hydro, new_dt_cooling, new_dt_grav);
  new_dt = min4(new_dt, new_dt_mhd, new_dt_chemistry, new_dt_forcing);

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
  integertime_t new_dti = make_integer_timestep(
      new_dt, p->time_bin, p->limiter_data.min_ngb_time_bin, e->ti_current,
      e->time_base_inv);

  if (e->policy & engine_policy_rt) {
    if (new_dti_rt <= new_dti) {
      /* enforce dt_hydro <= nsubcycles * dt_rt. The rare case where
       * new_dti_rt > new_dti will be handled in the parent function
       * that calls this one. */
      integertime_t max_subcycles = max(e->max_nr_rt_subcycles, 1);
      if (max_nr_timesteps / max_subcycles < new_dti_rt) {
        /* multiplication new_dti_rt * max_subcycles would overflow. This can
         * happen in rare cases, especially if the total physical time the
         * simulation should cover is small. So limit max_subcycles to a
         * reasonable value.
         * First find an integer guess for the maximal permissible number
         * of sub-cycles. Then find highest power-of-two below that guess.
         * Divide the guess by a factor of 2 to simplify the subsequent while
         * loop. The max() is there to prevent bad things happening. */
        const integertime_t max_subcycles_guess =
            max(1LL, max_nr_timesteps / (new_dti_rt * 2LL));
        max_subcycles = 1LL;
        while (max_subcycles_guess > max_subcycles) max_subcycles *= 2LL;
      }
      new_dti = min(new_dti, new_dti_rt * max_subcycles);
    }
  }

  return new_dti;
}

/**
 * @brief Compute the new (integer) time-step of a given #part
 *
 * @param p The #part.
 * @param xp The #xpart partner of p.
 * @param e The #engine (used to get some constants).
 */
__attribute__((always_inline)) INLINE static integertime_t get_part_rt_timestep(
    const struct part *restrict p, const struct xpart *restrict xp,
    const struct engine *restrict e) {

  if (!(e->policy & engine_policy_rt))
    return get_integer_timestep(num_time_bins);

  float new_dt =
      rt_compute_timestep(p, xp, e->rt_props, e->cosmology, e->hydro_properties,
                          e->physical_constants, e->internal_units);

  if ((e->policy & engine_policy_cosmology))
    /* Apply the maximal displacement constraint (FLT_MAX if non-cosmological)*/
    new_dt = min(new_dt, e->dt_max_RMS_displacement);

  /* Apply cosmology correction (This is 1 if non-cosmological) */
  new_dt *= e->cosmology->time_step_factor;

  /* Limit timestep within the allowed range */
  new_dt = min(new_dt, e->dt_max);

#ifdef SWIFT_RT_DEBUG_CHECKS
  /* Proper error will be caught in get_part_timestep(), so keep this as
   * debugging check only. */
  const float f = (float)max(e->max_nr_rt_subcycles, 1);
  if (new_dt < e->dt_min / f)
    error(
        "part (id=%lld) wants an RT time-step (%e) below dt_min/nr_subcycles "
        "(%e)",
        p->id, new_dt, e->dt_min / f);
#endif

  const integertime_t new_dti = make_integer_timestep(
      new_dt, p->rt_time_data.time_bin, p->rt_time_data.min_ngb_time_bin,
      e->ti_current, e->time_base_inv);

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
  float new_dt_stars = stars_compute_timestep(
      sp, e->stars_properties, (e->policy & engine_policy_cosmology),
      e->cosmology, e->time);

  /* Gravity time-step */
  float new_dt_self = FLT_MAX, new_dt_ext = FLT_MAX;

  if (e->policy & engine_policy_external_gravity)
    new_dt_ext = external_gravity_timestep(e->time, e->external_potential,
                                           e->physical_constants, sp->gpart);

  const float a_hydro[3] = {0.f, 0.f, 0.f};
  if (e->policy & engine_policy_self_gravity)
    new_dt_self = gravity_compute_timestep_self(
        sp->gpart, a_hydro, e->gravity_properties, e->cosmology);

  float new_dt_rt = FLT_MAX;
  if (e->policy & engine_policy_rt)
    new_dt_rt = rt_compute_spart_timestep(sp, e->rt_props, e->cosmology);

  /* Take the minimum of all */
  float new_dt = min4(new_dt_stars, new_dt_self, new_dt_ext, new_dt_rt);

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
      new_dt, sp->time_bin, num_time_bins, e->ti_current, e->time_base_inv);

  return new_dti;
}

/**
 * @brief Compute the new (integer) time-step of a given #bpart
 *
 * @param bp The #bpart.
 * @param e The #engine (used to get some constants).
 */
__attribute__((always_inline)) INLINE static integertime_t get_bpart_timestep(
    const struct bpart *restrict bp, const struct engine *restrict e) {

  /* Black hole internal time-step */
  float new_dt_black_holes = black_holes_compute_timestep(
      bp, e->black_holes_properties, e->physical_constants, e->cosmology);

  /* Gravity time-step */
  float new_dt_self = FLT_MAX, new_dt_ext = FLT_MAX;

  if (e->policy & engine_policy_external_gravity)
    new_dt_ext = external_gravity_timestep(e->time, e->external_potential,
                                           e->physical_constants, bp->gpart);

  const float a_hydro[3] = {0.f, 0.f, 0.f};
  if (e->policy & engine_policy_self_gravity)
    new_dt_self = gravity_compute_timestep_self(
        bp->gpart, a_hydro, e->gravity_properties, e->cosmology);

  /* Take the minimum of all */
  float new_dt = min3(new_dt_black_holes, new_dt_self, new_dt_ext);

  /* Apply the maximal dibslacement constraint (FLT_MAX  if non-cosmological)*/
  new_dt = min(new_dt, e->dt_max_RMS_displacement);

  /* Apply cosmology correction (This is 1 if non-cosmological) */
  new_dt *= e->cosmology->time_step_factor;

  /* Limit timestep within the allowed range */
  new_dt = min(new_dt, e->dt_max);
  if (new_dt < e->dt_min) {
    error("bpart (id=%lld) wants a time-step (%e) below dt_min (%e)", bp->id,
          new_dt, e->dt_min);
  }

  /* Convert to integer time */
  const integertime_t new_dti = make_integer_timestep(
      new_dt, bp->time_bin, num_time_bins, e->ti_current, e->time_base_inv);

  return new_dti;
}

/**
 * @brief Compute the new (integer) time-step of a given #sink.
 *
 * @param sink The #sink.
 * @param e The #engine (used to get some constants).
 */
__attribute__((always_inline)) INLINE static integertime_t get_sink_timestep(
    const struct sink *restrict sink, const struct engine *restrict e) {

  /* Sink time-step */
  float new_dt_sink = sink_compute_timestep(sink);

  /* Gravity time-step */
  float new_dt_self = FLT_MAX, new_dt_ext = FLT_MAX;

  if (e->policy & engine_policy_external_gravity)
    new_dt_ext = external_gravity_timestep(e->time, e->external_potential,
                                           e->physical_constants, sink->gpart);

  const float a_hydro[3] = {0.f, 0.f, 0.f};
  if (e->policy & engine_policy_self_gravity)
    new_dt_self = gravity_compute_timestep_self(
        sink->gpart, a_hydro, e->gravity_properties, e->cosmology);

  /* Take the minimum of all */
  float new_dt = min3(new_dt_sink, new_dt_self, new_dt_ext);

  /* Apply the maximal dibslacement constraint (FLT_MAX  if non-cosmological)*/
  new_dt = min(new_dt, e->dt_max_RMS_displacement);

  /* Apply cosmology correction (This is 1 if non-cosmological) */
  new_dt *= e->cosmology->time_step_factor;

  /* Limit timestep within the allowed range */
  new_dt = min(new_dt, e->dt_max);
  if (new_dt < e->dt_min) {
    error("sink (id=%lld) wants a time-step (%e) below dt_min (%e)", sink->id,
          new_dt, e->dt_min);
  }

  /* Convert to integer time */
  const integertime_t new_dti = make_integer_timestep(
      new_dt, sink->time_bin, num_time_bins, e->ti_current, e->time_base_inv);

  return new_dti;
}

#endif /* SWIFT_TIMESTEP_H */
