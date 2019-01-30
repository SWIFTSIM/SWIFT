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
#include "stars.h"

/**
 * @brief Perform the 'drift' operation on a #gpart.
 *
 * @param gp The #gpart to drift.
 * @param dt_drift The drift time-step.
 * @param ti_old Integer start of time-step (for debugging checks).
 * @param ti_current Integer end of time-step (for debugging checks).
 */
__attribute__((always_inline)) INLINE static void drift_gpart(
    struct gpart *restrict gp, double dt_drift, integertime_t ti_old,
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
  gp->x[0] += gp->v_full[0] * dt_drift;
  gp->x[1] += gp->v_full[1] * dt_drift;
  gp->x[2] += gp->v_full[2] * dt_drift;
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
 */
__attribute__((always_inline)) INLINE static void drift_part(
    struct part *restrict p, struct xpart *restrict xp, double dt_drift,
    double dt_kick_hydro, double dt_kick_grav, double dt_therm,
    integertime_t ti_old, integertime_t ti_current) {

#ifdef SWIFT_DEBUG_CHECKS
  if (p->ti_drift != ti_old)
    error(
        "particle has not been drifted to the current time p->ti_drift=%lld, "
        "c->ti_old=%lld, ti_current=%lld",
        p->ti_drift, ti_old, ti_current);

  p->ti_drift = ti_current;
#endif

  /* Drift... */
  p->x[0] += xp->v_full[0] * dt_drift;
  p->x[1] += xp->v_full[1] * dt_drift;
  p->x[2] += xp->v_full[2] * dt_drift;

  /* Predict velocities (for hydro terms) */
  p->v[0] += p->a_hydro[0] * dt_kick_hydro;
  p->v[1] += p->a_hydro[1] * dt_kick_hydro;
  p->v[2] += p->a_hydro[2] * dt_kick_hydro;

  p->v[0] += xp->a_grav[0] * dt_kick_grav;
  p->v[1] += xp->a_grav[1] * dt_kick_grav;
  p->v[2] += xp->a_grav[2] * dt_kick_grav;

  /* Predict the values of the extra fields */
  hydro_predict_extra(p, xp, dt_drift, dt_therm);

  /* Compute offsets since last cell construction */
  for (int k = 0; k < 3; k++) {
    const float dx = xp->v_full[k] * dt_drift;
    xp->x_diff[k] -= dx;
    xp->x_diff_sort[k] -= dx;
  }
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
}

#endif /* SWIFT_DRIFT_H */
