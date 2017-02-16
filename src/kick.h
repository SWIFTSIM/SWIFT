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
#include "stars.h"
#include "timeline.h"

/**
 * @brief Perform the 'kick' operation on a #gpart
 *
 * @param gp The #gpart to kick.
 * @param ti_start The starting (integer) time of the kick
 * @param ti_end The ending (integer) time of the kick
 * @param timeBase The minimal allowed time-step size.
 */
__attribute__((always_inline)) INLINE static void kick_gpart(
    struct gpart *restrict gp, integertime_t ti_start, integertime_t ti_end,
    double timeBase) {

  /* Time interval for this half-kick */
  const float dt = (ti_end - ti_start) * timeBase;

#ifdef SWIFT_DEBUG_CHECKS
  if (gp->ti_kick != ti_start)
    error(
        "g-particle has not been kicked to the current time gp->ti_kick=%lld, "
        "ti_start=%lld, ti_end=%lld",
        gp->ti_kick, ti_start, ti_end);

  gp->ti_kick = ti_end;
#endif

  /* Kick particles in momentum space */
  gp->v_full[0] += gp->a_grav[0] * dt;
  gp->v_full[1] += gp->a_grav[1] * dt;
  gp->v_full[2] += gp->a_grav[2] * dt;

  /* Kick extra variables */
  gravity_kick_extra(gp, dt);
}

/**
 * @brief Perform the 'kick' operation on a #part
 *
 * @param p The #part to kick.
 * @param xp The #xpart of the particle.
 * @param ti_start The starting (integer) time of the kick
 * @param ti_end The ending (integer) time of the kick
 * @param timeBase The minimal allowed time-step size.
 */
__attribute__((always_inline)) INLINE static void kick_part(
    struct part *restrict p, struct xpart *restrict xp, integertime_t ti_start,
    integertime_t ti_end, double timeBase) {

  /* Time interval for this half-kick */
  const float dt = (ti_end - ti_start) * timeBase;

#ifdef SWIFT_DEBUG_CHECKS
  if (p->ti_kick != ti_start)
    error(
        "particle has not been kicked to the current time p->ti_kick=%lld, "
        "ti_start=%lld, ti_end=%lld",
        p->ti_kick, ti_start, ti_end);

  p->ti_kick = ti_end;
#endif

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

  /* Extra kick work */
  hydro_kick_extra(p, xp, dt);
  if (p->gpart != NULL) gravity_kick_extra(p->gpart, dt);
}

/**
 * @brief Perform the 'kick' operation on a #spart
 *
 * @param sp The #spart to kick.
 * @param ti_start The starting (integer) time of the kick
 * @param ti_end The ending (integer) time of the kick
 * @param timeBase The minimal allowed time-step size.
 */
__attribute__((always_inline)) INLINE static void kick_spart(
    struct spart *restrict sp, integertime_t ti_start, integertime_t ti_end,
    double timeBase) {

  /* Time interval for this half-kick */
  const float dt = (ti_end - ti_start) * timeBase;

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->ti_kick != ti_start)
    error(
        "s-particle has not been kicked to the current time sp->ti_kick=%lld, "
        "ti_start=%lld, ti_end=%lld",
        sp->ti_kick, ti_start, ti_end);

  sp->ti_kick = ti_end;
#endif

  /* Acceleration from gravity */
  const float a[3] = {sp->gpart->a_grav[0], sp->gpart->a_grav[1],
                      sp->gpart->a_grav[2]};

  /* Kick particles in momentum space */
  sp->v[0] += a[0] * dt;
  sp->v[1] += a[1] * dt;
  sp->v[2] += a[2] * dt;
  sp->gpart->v_full[0] = sp->v[0];
  sp->gpart->v_full[1] = sp->v[1];
  sp->gpart->v_full[2] = sp->v[2];

  /* Kick extra variables */
  star_kick_extra(sp, dt);
}

#endif /* SWIFT_KICK_H */
