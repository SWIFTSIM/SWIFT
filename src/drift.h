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
#ifndef SWIFT_DRIFT_H
#define SWIFT_DRIFT_H

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "const.h"
#include "debug.h"
#include "dimension.h"
#include "hydro.h"
#include "part.h"

/**
 * @brief Perform the 'drift' operation on a #gpart
 *
 * @param gp The #gpart to drift.
 * @param dt The drift time-step
 * @param timeBase The minimal allowed time-step size.
 * @param ti_old Integer start of time-step
 * @param ti_current Integer end of time-step
 */
__attribute__((always_inline)) INLINE static void drift_gpart(
    struct gpart *restrict gp, double dt, double timeBase, integertime_t ti_old,
    integertime_t ti_current) {

#ifdef SWIFT_DEBUG_CHECKS
  if (gp->ti_drift != ti_old)
    error(
        "g-particle has not been drifted to the current time "
        "gp->ti_drift=%lld, "
        "c->ti_old=%lld, ti_current=%lld",
        gp->ti_drift, ti_old, ti_current);

  gp->ti_drift = ti_current;
#endif

  /* Drift... */
  gp->x[0] += gp->v_full[0] * dt;
  gp->x[1] += gp->v_full[1] * dt;
  gp->x[2] += gp->v_full[2] * dt;

  /* Compute offset since last cell construction */
  gp->x_diff[0] -= gp->v_full[0] * dt;
  gp->x_diff[1] -= gp->v_full[1] * dt;
  gp->x_diff[2] -= gp->v_full[2] * dt;
}

/**
 * @brief Perform the 'drift' operation on a #part
 *
 * @param p The #part to drift.
 * @param xp The #xpart of the particle.
 * @param dt The drift time-step
 * @param timeBase The minimal allowed time-step size.
 * @param ti_old Integer start of time-step
 * @param ti_current Integer end of time-step
 */
__attribute__((always_inline)) INLINE static void drift_part(
    struct part *restrict p, struct xpart *restrict xp, double dt,
    double timeBase, integertime_t ti_old, integertime_t ti_current) {

#ifdef SWIFT_DEBUG_CHECKS
  if (p->ti_drift != ti_old)
    error(
        "particle has not been drifted to the current time p->ti_drift=%lld, "
        "c->ti_old=%lld, ti_current=%lld",
        p->ti_drift, ti_old, ti_current);

  p->ti_drift = ti_current;
#endif

  /* Drift... */
  p->x[0] += xp->v_full[0] * dt;
  p->x[1] += xp->v_full[1] * dt;
  p->x[2] += xp->v_full[2] * dt;

  /* Predict velocities (for hydro terms) */
  p->v[0] += p->a_hydro[0] * dt;
  p->v[1] += p->a_hydro[1] * dt;
  p->v[2] += p->a_hydro[2] * dt;

  /* Predict the values of the extra fields */
  hydro_predict_extra(p, xp, dt);

  /* Compute offsets since last cell construction */
  for (int k = 0; k < 3; k++) {
    const float dx = xp->v_full[k] * dt;
    xp->x_diff[k] -= dx;
    xp->x_diff_sort[k] -= dx;
  }
}

/**
 * @brief Perform the 'drift' operation on a #spart
 *
 * @param sp The #spart to drift.
 * @param dt The drift time-step
 * @param timeBase The minimal allowed time-step size.
 * @param ti_old Integer start of time-step
 * @param ti_current Integer end of time-step
 */
__attribute__((always_inline)) INLINE static void drift_spart(
    struct spart *restrict sp, double dt, double timeBase, integertime_t ti_old,
    integertime_t ti_current) {

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->ti_drift != ti_old)
    error(
        "s-particle has not been drifted to the current time "
        "sp->ti_drift=%lld, "
        "c->ti_old=%lld, ti_current=%lld",
        sp->ti_drift, ti_old, ti_current);

  sp->ti_drift = ti_current;
#endif

  /* Drift... */
  sp->x[0] += sp->v[0] * dt;
  sp->x[1] += sp->v[1] * dt;
  sp->x[2] += sp->v[2] * dt;
}

#endif /* SWIFT_DRIFT_H */
