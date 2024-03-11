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
#ifndef SWIFT_DRIFT_H
#define SWIFT_DRIFT_H

/* Config parameters. */
#include <config.h>

/* Local headers. */
#include "adaptive_softening.h"
#include "black_holes.h"
#include "const.h"
#include "debug.h"
#include "dimension.h"
#include "engine.h"
#include "entropy_floor.h"
#include "hydro.h"
#include "hydro_properties.h"
#include "lightcone/lightcone_crossing.h"
#include "lightcone/lightcone_replications.h"
#include "mhd.h"
#include "part.h"
#include "rt.h"
#include "sink.h"
#include "stars.h"

/**
 * @brief Perform the 'drift' operation on a #gpart.
 *
 * @param gp The #gpart to drift.
 * @param dt_drift The drift time-step.
 * @param ti_old Integer start of time-step (for debugging checks).
 * @param ti_current Integer end of time-step (for debugging checks).
 * @param grav_props The properties of the gravity scheme.
 * @param e the #engine
 */
__attribute__((always_inline)) INLINE static void drift_gpart(
    struct gpart *restrict gp, double dt_drift, integertime_t ti_old,
    integertime_t ti_current, const struct gravity_props *grav_props,
    const struct engine *e, struct replication_list *replication_list,
    const double cell_loc[3]) {

#ifdef SWIFT_DEBUG_CHECKS
  if (gp->time_bin == time_bin_not_created) {
    error("Found an extra gpart in the drift");
  }

  if (gp->ti_drift != ti_old)
    error(
        "g-particle has not been drifted to the current time "
        "gp->ti_drift=%lld, "
        "c->ti_old=%lld, ti_current=%lld",
        gp->ti_drift, ti_old, ti_current);

  gp->ti_drift = ti_current;
#endif

#ifdef SWIFT_FIXED_BOUNDARY_PARTICLES

  /* Get the ID of the gpart */
  long long id = 0;
  if (gp->type == swift_type_gas)
    id = e->s->parts[-gp->id_or_neg_offset].id;
  else if (gp->type == swift_type_stars)
    id = e->s->sparts[-gp->id_or_neg_offset].id;
  else if (gp->type == swift_type_black_hole)
    id = e->s->bparts[-gp->id_or_neg_offset].id;
  else
    id = gp->id_or_neg_offset;

  /* Cancel the velocity of the particles */
  if (id < SWIFT_FIXED_BOUNDARY_PARTICLES) {

    /* Don't move! */
    gp->v_full[0] = 0.f;
    gp->v_full[1] = 0.f;
    gp->v_full[2] = 0.f;
  }
#endif

#ifdef WITH_LIGHTCONE
  /* Store initial position and velocity for lightcone check after the drift */
  const double x[3] = {gp->x[0], gp->x[1], gp->x[2]};
  const float v_full[3] = {gp->v_full[0], gp->v_full[1], gp->v_full[2]};
#endif

  /* Drift... */
  gp->x[0] += gp->v_full[0] * dt_drift;
  gp->x[1] += gp->v_full[1] * dt_drift;
  gp->x[2] += gp->v_full[2] * dt_drift;

  gravity_predict_extra(gp, grav_props);

#ifdef WITH_LIGHTCONE
  /* Check for lightcone crossing */
  switch (gp->type) {
    case swift_type_dark_matter:
    case swift_type_dark_matter_background:
    case swift_type_neutrino:
      /* This particle has no *part counterpart, so check for lightcone crossing
       * here */
      lightcone_check_particle_crosses(e, replication_list, x, v_full, gp,
                                       dt_drift, ti_old, ti_current, cell_loc);
      break;
    default:
      /* Particle has a counterpart or is of a type not supported in lightcones
       */
      break;
  }
#endif
}

/**
 * @brief Perform the 'drift' operation on a #part
 *
 * @param p The #part to drift.
 * @param xp The #xpart of the particle.
 * @param dt_drift The drift time-step
 * @param dt_kick_grav The kick time-step for gravity accelerations.
 * @param dt_kick_hydro The kick time-step for hydro accelerations.
 * @param dt_therm The drift time-step for thermodynamic quantities.
 * @param ti_old Integer start of time-step (for debugging checks).
 * @param ti_current Integer end of time-step (for debugging checks).
 * @param cosmo The cosmological model.
 * @param hydro_props The properties of the hydro scheme.
 * @param floor The properties of the entropy floor.
 */
__attribute__((always_inline)) INLINE static void drift_part(
    struct part *restrict p, struct xpart *restrict xp, double dt_drift,
    double dt_kick_hydro, double dt_kick_grav, double dt_therm,
    integertime_t ti_old, integertime_t ti_current, const struct engine *e,
    struct replication_list *replication_list, const double cell_loc[3]) {

  const struct cosmology *cosmo = e->cosmology;
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct entropy_floor_properties *entropy_floor = e->entropy_floor;
  const struct pressure_floor_props *pressure_floor = e->pressure_floor_props;

#ifdef SWIFT_DEBUG_CHECKS
  if (p->ti_drift != ti_old)
    error(
        "particle has not been drifted to the current time p->ti_drift=%lld, "
        "c->ti_old=%lld, ti_current=%lld",
        p->ti_drift, ti_old, ti_current);

  p->ti_drift = ti_current;
#endif

#ifdef SWIFT_FIXED_BOUNDARY_PARTICLES

  /* Get the ID of the gpart */
  const long long id = p->id;

  /* Cancel the velocity of the particles */
  if (id < SWIFT_FIXED_BOUNDARY_PARTICLES) {

    /* Don't move! */
    xp->v_full[0] = 0.f;
    xp->v_full[1] = 0.f;
    xp->v_full[2] = 0.f;
  }
#endif

#ifdef WITH_LIGHTCONE
  /* Store initial position and velocity for lightcone check after the drift */
  const double x[3] = {p->x[0], p->x[1], p->x[2]};
  const float v_full[3] = {xp->v_full[0], xp->v_full[1], xp->v_full[2]};
#endif

  /* Drift... */
  p->x[0] += xp->v_full[0] * dt_drift;
  p->x[1] += xp->v_full[1] * dt_drift;
  p->x[2] += xp->v_full[2] * dt_drift;

  /* Predict velocities (for hydro terms) */
  p->v[0] += p->a_hydro[0] * dt_kick_hydro;
  p->v[1] += p->a_hydro[1] * dt_kick_hydro;
  p->v[2] += p->a_hydro[2] * dt_kick_hydro;

  /* Predict velocities (for gravity terms) */
  if (p->gpart != NULL) {
    p->v[0] += xp->a_grav[0] * dt_kick_grav;
    p->v[1] += xp->a_grav[1] * dt_kick_grav;
    p->v[2] += xp->a_grav[2] * dt_kick_grav;
  }

  /* Predict the values of the extra fields */
  hydro_predict_extra(p, xp, dt_drift, dt_therm, dt_kick_grav, cosmo,
                      hydro_props, entropy_floor, pressure_floor);
  mhd_predict_extra(p, xp, dt_drift, dt_therm, cosmo, hydro_props,
                    entropy_floor);
  rt_predict_extra(p, xp, dt_drift);
  if (p->gpart) gravity_update_softening(p->gpart, p, e->gravity_properties);

  /* Compute offsets since last cell construction */
  for (int k = 0; k < 3; k++) {
    const float dx = xp->v_full[k] * dt_drift;
    xp->x_diff[k] -= dx;
    xp->x_diff_sort[k] -= dx;
  }

#ifdef SWIFT_FIXED_BOUNDARY_PARTICLES

  /* Cancel the velocity of the particles */
  if (id < SWIFT_FIXED_BOUNDARY_PARTICLES) {
    p->v[0] = 0.f;
    p->v[1] = 0.f;
    p->v[2] = 0.f;
  }
#endif

#ifdef WITH_LIGHTCONE
  /* Check if the particle crossed the lightcone */
  if (p->gpart)
    lightcone_check_particle_crosses(e, replication_list, x, v_full, p->gpart,
                                     dt_drift, ti_old, ti_current, cell_loc);
#endif
}

/**
 * @brief Perform the 'drift' operation on a #spart
 *
 * @param sp The #spart to drift.
 * @param dt_drift The drift time-step.
 * @param ti_old Integer start of time-step (for debugging checks).
 * @param ti_current Integer end of time-step (for debugging checks).
 */
__attribute__((always_inline)) INLINE static void drift_spart(
    struct spart *restrict sp, double dt_drift, integertime_t ti_old,
    integertime_t ti_current, const struct engine *e,
    struct replication_list *replication_list, const double cell_loc[3]) {

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->ti_drift != ti_old)
    error(
        "s-particle has not been drifted to the current time "
        "sp->ti_drift=%lld, "
        "c->ti_old=%lld, ti_current=%lld",
        sp->ti_drift, ti_old, ti_current);

  sp->ti_drift = ti_current;
#endif

#ifdef SWIFT_FIXED_BOUNDARY_PARTICLES

  /* Get the ID of the gpart */
  const long long id = sp->id;

  /* Cancel the velocity of the particles */
  if (id < SWIFT_FIXED_BOUNDARY_PARTICLES) {

    /* Don't move! */
    sp->v[0] = 0.f;
    sp->v[1] = 0.f;
    sp->v[2] = 0.f;
  }
#endif

#ifdef WITH_LIGHTCONE
  /* Store initial position and velocity for lightcone check after the drift */
  const double x[3] = {sp->x[0], sp->x[1], sp->x[2]};
  const float v_full[3] = {sp->v[0], sp->v[1], sp->v[2]};
#endif

  /* Drift... */
  sp->x[0] += sp->v[0] * dt_drift;
  sp->x[1] += sp->v[1] * dt_drift;
  sp->x[2] += sp->v[2] * dt_drift;

  /* Predict the values of the extra fields */
  stars_predict_extra(sp, dt_drift);

  /* Compute offsets since last cell construction */
  for (int k = 0; k < 3; k++) {
    const float dx = sp->v[k] * dt_drift;
    sp->x_diff[k] -= dx;
    sp->x_diff_sort[k] -= dx;
  }

#ifdef WITH_LIGHTCONE
  /* Check for lightcone crossing */
  if (sp->gpart)
    lightcone_check_particle_crosses(e, replication_list, x, v_full, sp->gpart,
                                     dt_drift, ti_old, ti_current, cell_loc);
#endif
}

/**
 * @brief Perform the 'drift' operation on a #bpart
 *
 * @param bp The #bpart to drift.
 * @param dt_drift The drift time-step.
 * @param ti_old Integer start of time-step (for debugging checks).
 * @param ti_current Integer end of time-step (for debugging checks).
 */
__attribute__((always_inline)) INLINE static void drift_bpart(
    struct bpart *restrict bp, double dt_drift, integertime_t ti_old,
    integertime_t ti_current, const struct engine *e,
    struct replication_list *replication_list, const double cell_loc[3]) {

#ifdef SWIFT_DEBUG_CHECKS
  if (bp->ti_drift != ti_old)
    error(
        "b-particle has not been drifted to the current time "
        "bp->ti_drift=%lld, "
        "c->ti_old=%lld, ti_current=%lld",
        bp->ti_drift, ti_old, ti_current);

  bp->ti_drift = ti_current;
#endif

#ifdef SWIFT_FIXED_BOUNDARY_PARTICLES

  /* Get the ID of the gpart */
  const long long id = bp->id;

  /* Cancel the velocity of the particles */
  if (id < SWIFT_FIXED_BOUNDARY_PARTICLES) {

    /* Don't move! */
    bp->v[0] = 0.f;
    bp->v[1] = 0.f;
    bp->v[2] = 0.f;
  }
#endif

#ifdef WITH_LIGHTCONE
  /* Store initial position and velocity for lightcone check after the drift */
  const double x[3] = {bp->x[0], bp->x[1], bp->x[2]};
  const float v_full[3] = {bp->v[0], bp->v[1], bp->v[2]};
#endif

  /* Drift... */
  bp->x[0] += bp->v[0] * dt_drift;
  bp->x[1] += bp->v[1] * dt_drift;
  bp->x[2] += bp->v[2] * dt_drift;

  /* Predict the values of the extra fields */
  black_holes_predict_extra(bp, dt_drift);

  /* Compute offsets since last cell construction */
  for (int k = 0; k < 3; k++) {
    const float dx = bp->v[k] * dt_drift;
    bp->x_diff[k] -= dx;
  }

#ifdef WITH_LIGHTCONE
  /* Check for lightcone crossing */
  if (bp->gpart)
    lightcone_check_particle_crosses(e, replication_list, x, v_full, bp->gpart,
                                     dt_drift, ti_old, ti_current, cell_loc);
#endif
}

/**
 * @brief Perform the 'drift' operation on a #sink
 *
 * @param sink The #sink to drift.
 * @param dt_drift The drift time-step.
 * @param ti_old Integer start of time-step (for debugging checks).
 * @param ti_current Integer end of time-step (for debugging checks).
 */
__attribute__((always_inline)) INLINE static void drift_sink(
    struct sink *restrict sink, double dt_drift, integertime_t ti_old,
    integertime_t ti_current) {

#ifdef SWIFT_DEBUG_CHECKS
  if (sink->ti_drift != ti_old) {
    error(
        "s-particle has not been drifted to the current time "
        "sink->ti_drift=%lld, "
        "c->ti_old=%lld, ti_current=%lld",
        sink->ti_drift, ti_old, ti_current);
  }

  sink->ti_drift = ti_current;
#endif

  /* Drift... */
  sink->x[0] += sink->v[0] * dt_drift;
  sink->x[1] += sink->v[1] * dt_drift;
  sink->x[2] += sink->v[2] * dt_drift;

  /* Predict the values of the extra fields */
  sink_predict_extra(sink, dt_drift);

  /* Compute offsets since last cell construction */
  for (int k = 0; k < 3; k++) {
    const float dx = sink->v[k] * dt_drift;
    sink->x_diff[k] -= dx;
  }
}

#endif /* SWIFT_DRIFT_H */
