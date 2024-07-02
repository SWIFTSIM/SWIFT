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
#include "black_holes.h"
#include "chemistry.h"
#include "engine.h"
#include "feedback.h"
#include "gravity.h"
#include "mhd.h"
#include "neutrino.h"
#include "particle_splitting.h"
#include "rt.h"
#include "sink.h"
#include "star_formation.h"
#include "stars.h"
#include "threadpool.h"
#include "tracers.h"

void space_first_init_parts_mapper(void *restrict map_data, int count,
                                   void *restrict extra_data) {

  struct part *restrict p = (struct part *)map_data;
  const struct space *restrict s = (struct space *)extra_data;
  const struct engine *e = s->e;

  const ptrdiff_t delta = p - s->parts;
  struct xpart *restrict xp = s->xparts + delta;

  /* Extract some constants */
  const struct cosmology *cosmo = s->e->cosmology;
  const struct phys_const *phys_const = s->e->physical_constants;
  const struct unit_system *us = s->e->internal_units;
  const float a_factor_vel = cosmo->a;

  const struct hydro_props *hydro_props = s->e->hydro_properties;
  const float u_init = hydro_props->initial_internal_energy;
  const float hydro_h_min_ratio = e->hydro_properties->h_min_ratio;

  const struct gravity_props *grav_props = s->e->gravity_properties;
  const int with_gravity = e->policy & engine_policy_self_gravity;

  const struct chemistry_global_data *chemistry = e->chemistry;
  const struct star_formation *star_formation = e->star_formation;
  const struct cooling_function_data *cool_func = e->cooling_func;
  const struct rt_props *rt_props = e->rt_props;

  /* Check that the smoothing lengths are non-zero */
  for (int k = 0; k < count; k++) {
    if (p[k].h <= 0.)
      error("Invalid value of smoothing length for part %lld h=%e", p[k].id,
            p[k].h);

    if (with_gravity) {
      const struct gpart *gp = p[k].gpart;
      const float softening = gravity_get_softening(gp, grav_props);
      p->h = max(p->h, softening * hydro_h_min_ratio);
    }
  }

  /* Convert velocities to internal units */
  for (int k = 0; k < count; k++) {
    p[k].v[0] *= a_factor_vel;
    p[k].v[1] *= a_factor_vel;
    p[k].v[2] *= a_factor_vel;

#ifdef HYDRO_DIMENSION_2D
    p[k].x[2] = 0.f;
    p[k].v[2] = 0.f;
#endif

#ifdef HYDRO_DIMENSION_1D
    p[k].x[1] = p[k].x[2] = 0.f;
    p[k].v[1] = p[k].v[2] = 0.f;
#endif
  }

  /* Overwrite the internal energy? */
  if (u_init > 0.f) {
    for (int k = 0; k < count; k++) {
      hydro_set_init_internal_energy(&p[k], u_init);
    }
  }

  /* Initialise the rest */
  for (int k = 0; k < count; k++) {

    hydro_first_init_part(&p[k], &xp[k]);
    mhd_first_init_part(&p[k], &xp[k], &hydro_props->mhd, s->dim[0]);
    p[k].limiter_data.min_ngb_time_bin = num_time_bins + 1;
    p[k].limiter_data.wakeup = time_bin_not_awake;
    p[k].limiter_data.to_be_synchronized = 0;

#ifdef WITH_CSDS
    csds_part_data_init(&xp[k].csds_data);
#endif

    /* Also initialise the chemistry */
    chemistry_first_init_part(phys_const, us, cosmo, chemistry, &p[k], &xp[k]);

    /* Also initialise the star formation */
    star_formation_first_init_part(phys_const, us, cosmo, star_formation, &p[k],
                                   &xp[k]);

    /* And the cooling */
    cooling_first_init_part(phys_const, us, hydro_props, cosmo, cool_func,
                            &p[k], &xp[k]);

    /* And the tracers */
    tracers_first_init_xpart(&p[k], &xp[k], us, phys_const, cosmo, hydro_props,
                             cool_func);

    /* And the black hole markers */
    black_holes_mark_part_as_not_swallowed(&p[k].black_holes_data);

    /* And the sink markers */
    sink_mark_part_as_not_swallowed(&p[k].sink_data);

    /* Also initialise the splitting data */
    particle_splitting_mark_part_as_not_split(&xp[k].split_data, p[k].id);

    /* And the radiative transfer */
    rt_first_init_timestep_data(&p[k]);
    rt_first_init_part(&p[k], cosmo, rt_props);

#ifdef SWIFT_DEBUG_CHECKS
    /* Check part->gpart->part linkeage. */
    if (p[k].gpart && p[k].gpart->id_or_neg_offset != -(k + delta))
      error("Invalid gpart -> part link");

    /* Initialise the time-integration check variables */
    p[k].ti_drift = 0;
    p[k].ti_kick = 0;
#endif
  }
}

/**
 * @brief Initialises all the particles by setting them into a valid state
 *
 * Calls hydro_first_init_part() on all the particles
 * Calls chemistry_first_init_part() on all the particles
 * Calls cooling_first_init_part() on all the particles
 */
void space_first_init_parts(struct space *s, int verbose) {

  const ticks tic = getticks();
  if (s->nr_parts > 0)
    threadpool_map(&s->e->threadpool, space_first_init_parts_mapper, s->parts,
                   s->nr_parts, sizeof(struct part), threadpool_auto_chunk_size,
                   s);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_first_init_gparts_mapper(void *restrict map_data, int count,
                                    void *restrict extra_data) {

  struct gpart *restrict gp = (struct gpart *)map_data;
  const struct space *restrict s = (struct space *)extra_data;

  const struct cosmology *cosmo = s->e->cosmology;
  const float a_factor_vel = cosmo->a;
  const struct gravity_props *grav_props = s->e->gravity_properties;

  /* Convert velocities to internal units */
  for (int k = 0; k < count; k++) {
    gp[k].v_full[0] *= a_factor_vel;
    gp[k].v_full[1] *= a_factor_vel;
    gp[k].v_full[2] *= a_factor_vel;

#ifdef HYDRO_DIMENSION_2D
    gp[k].x[2] = 0.f;
    gp[k].v_full[2] = 0.f;
#endif

#ifdef HYDRO_DIMENSION_1D
    gp[k].x[1] = gp[k].x[2] = 0.f;
    gp[k].v_full[1] = gp[k].v_full[2] = 0.f;
#endif
  }

  /* Initialise the rest */
  for (int k = 0; k < count; k++) {

    gravity_first_init_gpart(&gp[k], grav_props);

    if (gp[k].type == swift_type_neutrino)
      gravity_first_init_neutrino(&gp[k], s->e);

#ifdef WITH_CSDS
    csds_part_data_init(&gp[k].csds_data);
#endif

#ifdef SWIFT_DEBUG_CHECKS
    /* Initialise the time-integration check variables */
    gp[k].ti_drift = 0;
    gp[k].ti_kick = 0;
    gp[k].ti_kick_mesh = 0;
#endif
  }
}

/**
 * @brief Initialises all the g-particles by setting them into a valid state
 *
 * Calls gravity_first_init_gpart() on all the particles
 */
void space_first_init_gparts(struct space *s, int verbose) {

  const ticks tic = getticks();
  if (s->nr_gparts > 0)
    threadpool_map(&s->e->threadpool, space_first_init_gparts_mapper, s->gparts,
                   s->nr_gparts, sizeof(struct gpart),
                   threadpool_auto_chunk_size, s);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_first_init_sparts_mapper(void *restrict map_data, int count,
                                    void *restrict extra_data) {

  struct spart *restrict sp = (struct spart *)map_data;
  const struct space *restrict s = (struct space *)extra_data;
  const struct engine *e = s->e;

  const struct chemistry_global_data *chemistry = e->chemistry;

#ifdef SWIFT_DEBUG_CHECKS
  const ptrdiff_t delta = sp - s->sparts;
#endif

  const float initial_h = s->initial_spart_h;

  const int with_feedback = (e->policy & engine_policy_feedback);
  const int with_cosmology = (e->policy & engine_policy_cosmology);

  const struct cosmology *cosmo = e->cosmology;
  const struct stars_props *stars_properties = e->stars_properties;
  const struct feedback_props *feedback_properties = e->feedback_props;
  const struct phys_const *phys_const = s->e->physical_constants;
  const struct unit_system *us = s->e->internal_units;
  const float a_factor_vel = cosmo->a;

  /* Convert velocities to internal units */
  for (int k = 0; k < count; k++) {

    sp[k].v[0] *= a_factor_vel;
    sp[k].v[1] *= a_factor_vel;
    sp[k].v[2] *= a_factor_vel;

    /* Imposed smoothing length from parameter file */
    if (initial_h != -1.f) {
      sp[k].h = initial_h;
    }

#ifdef HYDRO_DIMENSION_2D
    sp[k].x[2] = 0.f;
    sp[k].v[2] = 0.f;
#endif

#ifdef HYDRO_DIMENSION_1D
    sp[k].x[1] = sp[k].x[2] = 0.f;
    sp[k].v[1] = sp[k].v[2] = 0.f;
#endif
  }

  /* Check that the smoothing lengths are non-zero */
  for (int k = 0; k < count; k++) {
    if (with_feedback && sp[k].h <= 0.)
      error("Invalid value of smoothing length for spart %lld h=%e", sp[k].id,
            sp[k].h);
  }

  /* Initialise the rest */
  for (int k = 0; k < count; k++) {

    stars_first_init_spart(&sp[k], stars_properties, with_cosmology, cosmo->a,
                           e->time);

#ifdef WITH_CSDS
    csds_part_data_init(&sp[k].csds_data);
#endif

    /* And the tracers */
    tracers_first_init_spart(&sp[k], us, phys_const, cosmo);

    /* Also initialise the chemistry */
    chemistry_first_init_spart(chemistry, &sp[k]);

    /* Also initialise the feedback */
    feedback_first_init_spart(&sp[k], feedback_properties);

    /* Also initialise the splitting data */
    particle_splitting_mark_part_as_not_split(&sp[k].split_data, sp[k].id);

    /* And radiative transfer data */
    rt_first_init_spart(&sp[k]);

#ifdef SWIFT_DEBUG_CHECKS
    if (sp[k].gpart && sp[k].gpart->id_or_neg_offset != -(k + delta))
      error("Invalid gpart -> spart link");

    /* Initialise the time-integration check variables */
    sp[k].ti_drift = 0;
    sp[k].ti_kick = 0;
#endif
  }
}

/**
 * @brief Initialises all the s-particles by setting them into a valid state
 *
 * Calls stars_first_init_spart() on all the particles
 */
void space_first_init_sparts(struct space *s, int verbose) {
  const ticks tic = getticks();
  if (s->nr_sparts > 0)
    threadpool_map(&s->e->threadpool, space_first_init_sparts_mapper, s->sparts,
                   s->nr_sparts, sizeof(struct spart),
                   threadpool_auto_chunk_size, s);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_first_init_bparts_mapper(void *restrict map_data, int count,
                                    void *restrict extra_data) {

  struct bpart *restrict bp = (struct bpart *)map_data;
  const struct space *restrict s = (struct space *)extra_data;
  const struct engine *e = s->e;
  const struct black_holes_props *props = e->black_holes_properties;

#ifdef SWIFT_DEBUG_CHECKS
  const ptrdiff_t delta = bp - s->bparts;
#endif

  const float initial_h = s->initial_bpart_h;

  const struct cosmology *cosmo = e->cosmology;
  const struct phys_const *phys_const = s->e->physical_constants;
  const struct unit_system *us = s->e->internal_units;
  const float a_factor_vel = cosmo->a;

  /* Convert velocities to internal units */
  for (int k = 0; k < count; k++) {

    bp[k].v[0] *= a_factor_vel;
    bp[k].v[1] *= a_factor_vel;
    bp[k].v[2] *= a_factor_vel;

    /* Imposed smoothing length from parameter file */
    if (initial_h != -1.f) {
      bp[k].h = initial_h;
    }

#ifdef HYDRO_DIMENSION_2D
    bp[k].x[2] = 0.f;
    bp[k].v[2] = 0.f;
#endif

#ifdef HYDRO_DIMENSION_1D
    bp[k].x[1] = bp[k].x[2] = 0.f;
    bp[k].v[1] = bp[k].v[2] = 0.f;
#endif
  }

  /* Check that the smoothing lengths are non-zero */
  for (int k = 0; k < count; k++) {
    if (bp[k].h <= 0.)
      error("Invalid value of smoothing length for bpart %lld h=%e", bp[k].id,
            bp[k].h);
  }

  /* Initialise the rest */
  for (int k = 0; k < count; k++) {

    black_holes_first_init_bpart(&bp[k], props);

    /* And the tracers */
    tracers_first_init_bpart(&bp[k], us, phys_const, cosmo);

    /* And the splitting data */
    particle_splitting_mark_part_as_not_split(&bp[k].split_data, bp[k].id);

    /* And the black hole merger markers */
    black_holes_mark_bpart_as_not_swallowed(&bp[k].merger_data);

#ifdef SWIFT_DEBUG_CHECKS
    if (bp[k].gpart && bp[k].gpart->id_or_neg_offset != -(k + delta))
      error("Invalid gpart -> bpart link");

    /* Initialise the time-integration check variables */
    bp[k].ti_drift = 0;
    bp[k].ti_kick = 0;
#endif
  }
}

/**
 * @brief Initialises all the b-particles by setting them into a valid state
 *
 * Calls stars_first_init_bpart() on all the particles
 */
void space_first_init_bparts(struct space *s, int verbose) {
  const ticks tic = getticks();
  if (s->nr_bparts > 0)
    threadpool_map(&s->e->threadpool, space_first_init_bparts_mapper, s->bparts,
                   s->nr_bparts, sizeof(struct bpart),
                   threadpool_auto_chunk_size, s);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}

void space_first_init_sinks_mapper(void *restrict map_data, int count,
                                   void *restrict extra_data) {

  struct sink *restrict sink = (struct sink *)map_data;
  const struct space *restrict s = (struct space *)extra_data;
  const struct engine *e = s->e;
  const struct sink_props *props = e->sink_properties;

#ifdef SWIFT_DEBUG_CHECKS
  const ptrdiff_t delta = sink - s->sinks;
#endif

  const struct cosmology *cosmo = e->cosmology;
  const float a_factor_vel = cosmo->a;

  /* Convert velocities to internal units */
  for (int k = 0; k < count; k++) {

    sink[k].v[0] *= a_factor_vel;
    sink[k].v[1] *= a_factor_vel;
    sink[k].v[2] *= a_factor_vel;

#ifdef HYDRO_DIMENSION_2D
    sink[k].x[2] = 0.f;
    sink[k].v[2] = 0.f;
#endif

#ifdef HYDRO_DIMENSION_1D
    sink[k].x[1] = sink[k].x[2] = 0.f;
    sink[k].v[1] = sink[k].v[2] = 0.f;
#endif
  }

  /* Initialise the rest */
  for (int k = 0; k < count; k++) {

    sink_first_init_sink(&sink[k], props, e);

#ifdef SWIFT_DEBUG_CHECKS
    if (sink[k].gpart && sink[k].gpart->id_or_neg_offset != -(k + delta))
      error("Invalid gpart -> sink link");

    /* Initialise the time-integration check variables */
    sink[k].ti_drift = 0;
    sink[k].ti_kick = 0;
#endif
  }
}

/**
 * @brief Initialises all the sink-particles by setting them into a valid state
 *
 * Calls stars_first_init_sink() on all the particles
 */
void space_first_init_sinks(struct space *s, int verbose) {
  const ticks tic = getticks();
  if (s->nr_sinks > 0)
    threadpool_map(&s->e->threadpool, space_first_init_sinks_mapper, s->sinks,
                   s->nr_sinks, sizeof(struct sink), threadpool_auto_chunk_size,
                   s);

  if (verbose)
    message("took %.3f %s.", clocks_from_ticks(getticks() - tic),
            clocks_getunit());
}
