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
#ifndef SWIFT_KICK_H
#define SWIFT_KICK_H

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "const.h"
#include "debug.h"
#include "timeline.h"

/**
 * @brief Perform the 'kick' operation on a #gpart
 *
 * @param gp The #gpart to kick.
 * @param new_dti The (integer) time-step for this kick.
 * @param timeBase The minimal allowed time-step size.
 */
__attribute__((always_inline)) INLINE static void kick_gpart(
    struct gpart *restrict gp, int new_dti, double timeBase) {

  /* Compute the time step for this kick */
  const int ti_start = (gp->ti_begin + gp->ti_end) / 2;
  const int ti_end = gp->ti_end + new_dti / 2;
  const float dt = (ti_end - ti_start) * timeBase;
  const float half_dt = (ti_end - gp->ti_end) * timeBase;

  /* Move particle forward in time */
  gp->ti_begin = gp->ti_end;
  gp->ti_end = gp->ti_begin + new_dti;

  /* Kick particles in momentum space */
  gp->v_full[0] += gp->a_grav[0] * dt;
  gp->v_full[1] += gp->a_grav[1] * dt;
  gp->v_full[2] += gp->a_grav[2] * dt;

  /* Extra kick work */
  gravity_kick_extra(gp, dt, half_dt);
}

/**
 * @brief Perform the 'kick' operation on a #part
 *
 * @param p The #part to kick.
 * @param xp The #xpart of the particle.
 * @param new_dti The (integer) time-step for this kick.
 * @param timeBase The minimal allowed time-step size.
 */
__attribute__((always_inline)) INLINE static void kick_part(
    struct part *restrict p, struct xpart *restrict xp, int new_dti,
    double timeBase) {

  /* Compute the time step for this kick */
  const integertime_t ti_begin = get_integer_time_begin(t1, p->time_bin);
  const integertime_t ti_end = get_integer_time_end(t1, p->time_bin);
  const int ti_start = (ti_begin + ti_end) / 2;
  const int ti_end = ti_end + new_dti / 2;
  const float dt = (ti_end - ti_start) * timeBase;
  const float half_dt = (ti_end - p->ti_end) * timeBase;

  /* Move particle forward in time */
  p->ti_begin = p->ti_end;
  p->ti_end = p->ti_begin + new_dti;
  if (p->gpart != NULL) {
    p->gpart->ti_begin = p->ti_begin;
    p->gpart->ti_end = p->ti_end;
  }

  /* Get the acceleration */
  float a_tot[3] = {p->a_hydro[0], p->a_hydro[1], p->a_hydro[2]};
  if (p->gpart != NULL) {
    a_tot[0] += p->gpart->a_grav[0];
    a_tot[1] += p->gpart->a_grav[1];
    a_tot[2] += p->gpart->a_grav[2];
  }

  /* Kick particles in momentum space */
  xp->v_full[0] += a_tot[0] * dt;
  xp->v_full[1] += a_tot[1] * dt;
  xp->v_full[2] += a_tot[2] * dt;
  if (p->gpart != NULL) {
    p->gpart->v_full[0] = xp->v_full[0];
    p->gpart->v_full[1] = xp->v_full[1];
    p->gpart->v_full[2] = xp->v_full[2];
  }

  /* Go back by half-step for the hydro velocity */
  p->v[0] = xp->v_full[0] - half_dt * a_tot[0];
  p->v[1] = xp->v_full[1] - half_dt * a_tot[1];
  p->v[2] = xp->v_full[2] - half_dt * a_tot[2];

  /* Extra kick work */
  hydro_kick_extra(p, xp, dt, half_dt);
  if (p->gpart != NULL) gravity_kick_extra(p->gpart, dt, half_dt);
}

#endif /* SWIFT_KICK_H */
