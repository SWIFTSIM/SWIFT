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
#ifndef SWIFT_KICK_H
#define SWIFT_KICK_H

/* Config parameters. */
#include <config.h>

/* Local headers. */
#include "black_holes.h"
#include "chemistry_additions.h"
#include "const.h"
#include "debug.h"
#include "mhd.h"
#include "rt.h"
#include "sink.h"
#include "stars.h"
#include "timeline.h"

/**
 * @brief Compute the time-step length for a gravity kick.
 *
 * @param ti_beg Start of the time-step (on the integer time-line).
 * @param ti_end End of the time-step (on the integer time-line).
 * @param time_base Minimal time-step size on the time-line.
 * @param with_cosmology Are we running with cosmology integration?
 * @param cosmo The #cosmology object.
 *
 * @return The time-step size for the gravity kick (internal units).
 */
__attribute__((always_inline)) INLINE static double kick_get_grav_kick_dt(
    const integertime_t ti_beg, const integertime_t ti_end,
    const double time_base, const int with_cosmology,
    const struct cosmology *cosmo) {

  if (with_cosmology) {
    return cosmology_get_grav_kick_factor(cosmo, ti_beg, ti_end);
  } else {
    return (ti_end - ti_beg) * time_base;
  }
}

/**
 * @brief Compute the time-step length for a hydro kick.
 *
 * @param ti_beg Start of the time-step (on the integer time-line).
 * @param ti_end End of the time-step (on the integer time-line).
 * @param time_base Minimal time-step size on the time-line.
 * @param with_cosmology Are we running with cosmology integration?
 * @param cosmo The #cosmology object.
 *
 * @return The time-step size for the hydro kick (internal units).
 */
__attribute__((always_inline)) INLINE static double kick_get_hydro_kick_dt(
    const integertime_t ti_beg, const integertime_t ti_end,
    const double time_base, const int with_cosmology,
    const struct cosmology *cosmo) {

  if (with_cosmology) {
    return cosmology_get_hydro_kick_factor(cosmo, ti_beg, ti_end);
  } else {
    return (ti_end - ti_beg) * time_base;
  }
}

/**
 * @brief Compute the time-step length for a thermal kick.
 *
 * @param ti_beg Start of the time-step (on the integer time-line).
 * @param ti_end End of the time-step (on the integer time-line).
 * @param time_base Minimal time-step size on the time-line.
 * @param with_cosmology Are we running with cosmology integration?
 * @param cosmo The #cosmology object.
 *
 * @return The time-step size for the thermal kick (internal units).
 */
__attribute__((always_inline)) INLINE static double kick_get_therm_kick_dt(
    const integertime_t ti_beg, const integertime_t ti_end,
    const double time_base, const int with_cosmology,
    const struct cosmology *cosmo) {

  if (with_cosmology) {
    return cosmology_get_therm_kick_factor(cosmo, ti_beg, ti_end);
  } else {
    return (ti_end - ti_beg) * time_base;
  }
}

/**
 * @brief Compute the time-step length for a gravity correction kick.
 *
 * @param ti_beg Start of the time-step (on the integer time-line).
 * @param ti_end End of the time-step (on the integer time-line).
 * @param time_base Minimal time-step size on the time-line.
 * @param with_cosmology Are we running with cosmology integration?
 * @param cosmo The #cosmology object.
 *
 * @return The time-step size for the gravity correction kick (internal units).
 */
__attribute__((always_inline)) INLINE static double kick_get_corr_kick_dt(
    const integertime_t ti_beg, const integertime_t ti_end,
    const double time_base, const int with_cosmology,
    const struct cosmology *cosmo) {

  if (with_cosmology) {
    return cosmology_get_corr_kick_factor(cosmo, ti_beg, ti_end);
  } else {
    return (ti_end - ti_beg) * time_base;
  }
}

/**
 * @brief Perform the 'kick' operation on a #gpart
 *
 * @param gp The #gpart to kick.
 * @param dt_kick_grav The kick time-step for gravity accelerations.
 * @param ti_start The starting (integer) time of the kick (for debugging
 * checks).
 * @param ti_end The ending (integer) time of the kick (for debugging checks).
 * @param dt_kick_mesh_grav The kick time-step for mesh gravity accelerations.
 * @param ti_start_mesh The starting (integer) time of the mesh kick (for
 * debugging checks).
 * @param ti_end_mesh The ending (integer) time of the mesh kick (for debugging
 * checks).
 */
__attribute__((always_inline)) INLINE static void kick_gpart(
    struct gpart *restrict gp, const double dt_kick_grav,
    const integertime_t ti_start, const integertime_t ti_end,
    const double dt_kick_mesh_grav, const integertime_t ti_start_mesh,
    const integertime_t ti_end_mesh) {

#ifdef SWIFT_DEBUG_CHECKS
  if (gp->time_bin == time_bin_not_created) {
    error("Found an extra gpart in the kick");
  }

  if (gp->ti_kick != ti_start)
    error(
        "g-particle has not been kicked to the current time gp->ti_kick=%lld, "
        "ti_start=%lld, ti_end=%lld id=%lld",
        gp->ti_kick, ti_start, ti_end, gp->id_or_neg_offset);

  gp->ti_kick = ti_end;

  if (ti_start_mesh != -1 && gp->ti_kick_mesh != ti_start_mesh)
    error(
        "g-particle has not been kicked (mesh) to the current time "
        "gp->ti_kick_mesh=%lld, "
        "ti_start_mesh=%lld, ti_end_mesh=%lld id=%lld",
        gp->ti_kick_mesh, ti_start_mesh, ti_end_mesh, gp->id_or_neg_offset);

  /* Record the mesh kick if we are doing one */
  if (ti_start_mesh != -1) gp->ti_kick_mesh = ti_end_mesh;

  if (ti_start_mesh == -1 && dt_kick_mesh_grav != 0.)
    error("Incorrect dt_kick for the mesh! %e (should be 0)",
          dt_kick_mesh_grav);

  if (ti_start_mesh != -1 && dt_kick_mesh_grav == 0.)
    error("Incorrect dt_kick for the mesh! %e (should not be 0)",
          dt_kick_mesh_grav);
#endif

  /* Kick particles in momentum space */
  gp->v_full[0] += gp->a_grav[0] * dt_kick_grav;
  gp->v_full[1] += gp->a_grav[1] * dt_kick_grav;
  gp->v_full[2] += gp->a_grav[2] * dt_kick_grav;

  /* Same for the long-range forces */
  gp->v_full[0] += gp->a_grav_mesh[0] * dt_kick_mesh_grav;
  gp->v_full[1] += gp->a_grav_mesh[1] * dt_kick_mesh_grav;
  gp->v_full[2] += gp->a_grav_mesh[2] * dt_kick_mesh_grav;

  /* Kick extra variables */
  gravity_kick_extra(gp, dt_kick_grav);
}

/**
 * @brief Perform the 'kick' operation on a #part
 *
 * @param p The #part to kick.
 * @param xp The #xpart of the particle.
 * @param dt_kick_hydro The kick time-step for hydro accelerations.
 * @param dt_kick_grav The kick time-step for gravity accelerations.
 * @param dt_kick_mesh_grav The kick time-step for mesh gravity accelerations.
 * @param dt_kick_therm The kick time-step for changes in thermal state.
 * @param dt_kick_corr The kick time-step for the gizmo-mfv gravity correction.
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme.
 * @param floor_props Properties of the entropy floor.
 * @param ti_start The starting (integer) time of the kick (for debugging
 * checks).
 * @param ti_end The ending (integer) time of the kick (for debugging checks).
 * @param ti_start_mesh The starting (integer) time of the mesh kick (for
 * debugging checks).
 * @param ti_end_mesh The ending (integer) time of the mesh kick (for debugging
 * checks).
 */
__attribute__((always_inline)) INLINE static void kick_part(
    struct part *restrict p, struct xpart *restrict xp,
    const double dt_kick_hydro, const double dt_kick_grav,
    const double dt_kick_mesh_grav, const double dt_kick_therm,
    const double dt_kick_corr, const struct cosmology *cosmo,
    const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props,
    const integertime_t ti_start, const integertime_t ti_end,
    const integertime_t ti_start_mesh, const integertime_t ti_end_mesh) {

#ifdef SWIFT_DEBUG_CHECKS
  if (p->ti_kick != ti_start)
    error(
        "particle has not been kicked to the current time p->ti_kick=%lld, "
        "ti_start=%lld, ti_end=%lld id=%lld time_bin=%d wakeup=%d",
        p->ti_kick, ti_start, ti_end, p->id, p->time_bin,
        p->limiter_data.wakeup);

  p->ti_kick = ti_end;

  if (ti_start_mesh == -1 && dt_kick_mesh_grav != 0.)
    error("Incorrect dt_kick for the mesh! %e (should be 0)",
          dt_kick_mesh_grav);

  if (ti_start_mesh != -1 && dt_kick_mesh_grav == 0.)
    error("Incorrect dt_kick for the mesh! %e (should not be 0)",
          dt_kick_mesh_grav);
#endif

  /* Kick particles in momentum space (hydro acc.) */
  xp->v_full[0] += p->a_hydro[0] * dt_kick_hydro;
  xp->v_full[1] += p->a_hydro[1] * dt_kick_hydro;
  xp->v_full[2] += p->a_hydro[2] * dt_kick_hydro;

  /* Kick particles in momentum space (grav acc.) */
  if (p->gpart != NULL) {
    xp->v_full[0] += p->gpart->a_grav[0] * dt_kick_grav;
    xp->v_full[1] += p->gpart->a_grav[1] * dt_kick_grav;
    xp->v_full[2] += p->gpart->a_grav[2] * dt_kick_grav;
  }

  /* Kick particles in momentum space (mesh grav acc.) */
  if (p->gpart != NULL) {
    xp->v_full[0] += p->gpart->a_grav_mesh[0] * dt_kick_mesh_grav;
    xp->v_full[1] += p->gpart->a_grav_mesh[1] * dt_kick_mesh_grav;
    xp->v_full[2] += p->gpart->a_grav_mesh[2] * dt_kick_mesh_grav;
  }

  /* Give the gpart friend the same velocity */
  if (p->gpart != NULL) {
    p->gpart->v_full[0] = xp->v_full[0];
    p->gpart->v_full[1] = xp->v_full[1];
    p->gpart->v_full[2] = xp->v_full[2];
  }

  /* Extra kick work (thermal quantities etc.) */
  /* for the GEAR RT, we need to do this before we update
   * the particle masses in hydro_kick_extra */
  rt_kick_extra(p, dt_kick_therm, dt_kick_grav, dt_kick_hydro, dt_kick_corr,
                cosmo, hydro_props);
  /* Similarly, we must apply the chemistry metal fluxes before updating the
   * particle masses */
  chemistry_kick_extra(p, dt_kick_therm, dt_kick_grav, dt_kick_hydro,
                       dt_kick_corr, cosmo, hydro_props);
  hydro_kick_extra(p, xp, dt_kick_therm, dt_kick_grav, dt_kick_mesh_grav,
                   dt_kick_hydro, dt_kick_corr, cosmo, hydro_props,
                   floor_props);
  mhd_kick_extra(p, xp, dt_kick_therm, dt_kick_grav, dt_kick_hydro,
                 dt_kick_corr, cosmo, hydro_props, floor_props);
  if (p->gpart != NULL) gravity_kick_extra(p->gpart, dt_kick_grav);
}

/**
 * @brief Perform the 'kick' operation on a #spart
 *
 * @param sp The #spart to kick.
 * @param dt_kick_grav The kick time-step for gravity accelerations.
 * @param ti_start The starting (integer) time of the kick (for debugging
 * checks).
 * @param ti_end The ending (integer) time of the kick (for debugging checks).
 * @param dt_kick_mesh_grav The kick time-step for mesh gravity accelerations.
 * @param ti_start_mesh The starting (integer) time of the mesh kick (for
 * debugging checks).
 * @param ti_end_mesh The ending (integer) time of the mesh kick (for debugging
 * checks).
 */
__attribute__((always_inline)) INLINE static void kick_spart(
    struct spart *restrict sp, const double dt_kick_grav,
    const integertime_t ti_start, const integertime_t ti_end,
    const double dt_kick_mesh_grav, const integertime_t ti_start_mesh,
    const integertime_t ti_end_mesh) {

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->ti_kick != ti_start)
    error(
        "s-particle has not been kicked to the current time sp->ti_kick=%lld, "
        "ti_start=%lld, ti_end=%lld id=%lld",
        sp->ti_kick, ti_start, ti_end, sp->id);

  sp->ti_kick = ti_end;

  if (ti_start_mesh == -1 && dt_kick_mesh_grav != 0.)
    error("Incorrect dt_kick for the mesh! %e (should be 0)",
          dt_kick_mesh_grav);

  if (ti_start_mesh != -1 && dt_kick_mesh_grav == 0.)
    error("Incorrect dt_kick for the mesh! %e (should not be 0)",
          dt_kick_mesh_grav);
#endif

  /* Kick particles in momentum space */
  sp->v[0] += sp->gpart->a_grav[0] * dt_kick_grav;
  sp->v[1] += sp->gpart->a_grav[1] * dt_kick_grav;
  sp->v[2] += sp->gpart->a_grav[2] * dt_kick_grav;

  /* Kick particles in momentum space (mesh forces) */
  sp->v[0] += sp->gpart->a_grav_mesh[0] * dt_kick_mesh_grav;
  sp->v[1] += sp->gpart->a_grav_mesh[1] * dt_kick_mesh_grav;
  sp->v[2] += sp->gpart->a_grav_mesh[2] * dt_kick_mesh_grav;

  /* Give the gpart friend the same velocity */
  sp->gpart->v_full[0] = sp->v[0];
  sp->gpart->v_full[1] = sp->v[1];
  sp->gpart->v_full[2] = sp->v[2];

  /* Kick extra variables */
  stars_kick_extra(sp, dt_kick_grav);
}

/**
 * @brief Perform the 'kick' operation on a #bpart
 *
 * @param bp The #bpart to kick.
 * @param dt_kick_grav The kick time-step for gravity accelerations.
 * @param ti_start The starting (integer) time of the kick (for debugging
 * checks).
 * @param ti_end The ending (integer) time of the kick (for debugging checks).
 * @param dt_kick_mesh_grav The kick time-step for mesh gravity accelerations.
 * @param ti_start_mesh The starting (integer) time of the mesh kick (for
 * debugging checks).
 * @param ti_end_mesh The ending (integer) time of the mesh kick (for debugging
 * checks).
 */
__attribute__((always_inline)) INLINE static void kick_bpart(
    struct bpart *restrict bp, const double dt_kick_grav,
    const integertime_t ti_start, const integertime_t ti_end,
    const double dt_kick_mesh_grav, const integertime_t ti_start_mesh,
    const integertime_t ti_end_mesh) {

#ifdef SWIFT_DEBUG_CHECKS
  if (bp->ti_kick != ti_start)
    error(
        "s-particle has not been kicked to the current time bp->ti_kick=%lld, "
        "ti_start=%lld, ti_end=%lld id=%lld",
        bp->ti_kick, ti_start, ti_end, bp->id);

  bp->ti_kick = ti_end;

  if (ti_start_mesh == -1 && dt_kick_mesh_grav != 0.)
    error("Incorrect dt_kick for the mesh! %e (should be 0)",
          dt_kick_mesh_grav);

  if (ti_start_mesh != -1 && dt_kick_mesh_grav == 0.)
    error("Incorrect dt_kick for the mesh! %e (should not be 0)",
          dt_kick_mesh_grav);
#endif

  /* Kick particles in momentum space */
  bp->v[0] += bp->gpart->a_grav[0] * dt_kick_grav;
  bp->v[1] += bp->gpart->a_grav[1] * dt_kick_grav;
  bp->v[2] += bp->gpart->a_grav[2] * dt_kick_grav;

  /* Kick particles in momentum space (mesh forces)*/
  bp->v[0] += bp->gpart->a_grav_mesh[0] * dt_kick_mesh_grav;
  bp->v[1] += bp->gpart->a_grav_mesh[1] * dt_kick_mesh_grav;
  bp->v[2] += bp->gpart->a_grav_mesh[2] * dt_kick_mesh_grav;

  /* Give the gpart friend the same velocity */
  bp->gpart->v_full[0] = bp->v[0];
  bp->gpart->v_full[1] = bp->v[1];
  bp->gpart->v_full[2] = bp->v[2];

  /* Kick extra variables */
  black_holes_kick_extra(bp, dt_kick_grav);
}

/**
 * @brief Perform the 'kick' operation on a #sink
 *
 * @param sink The #sink to kick.
 * @param dt_kick_grav The kick time-step for gravity accelerations.
 * @param ti_start The starting (integer) time of the kick (for debugging
 * checks).
 * @param ti_end The ending (integer) time of the kick (for debugging checks).
 * @param dt_kick_mesh_grav The kick time-step for mesh gravity accelerations.
 * @param ti_start_mesh The starting (integer) time of the mesh kick (for
 * debugging checks).
 * @param ti_end_mesh The ending (integer) time of the mesh kick (for debugging
 * checks).
 */
__attribute__((always_inline)) INLINE static void kick_sink(
    struct sink *restrict sink, const double dt_kick_grav,
    const integertime_t ti_start, const integertime_t ti_end,
    const double dt_kick_mesh_grav, const integertime_t ti_start_mesh,
    const integertime_t ti_end_mesh) {

#ifdef SWIFT_DEBUG_CHECKS
  if (sink->ti_kick != ti_start)
    error(
        "sink-particle has not been kicked to the current time "
        "sink->ti_kick=%lld, "
        "ti_start=%lld, ti_end=%lld id=%lld",
        sink->ti_kick, ti_start, ti_end, sink->id);

  sink->ti_kick = ti_end;

  if (ti_start_mesh == -1 && dt_kick_mesh_grav != 0.)
    error("Incorrect dt_kick for the mesh! %e (should be 0)",
          dt_kick_mesh_grav);

  if (ti_start_mesh != -1 && dt_kick_mesh_grav == 0.)
    error("Incorrect dt_kick for the mesh! %e (should not be 0)",
          dt_kick_mesh_grav);

#endif

  /* Kick particles in momentum space */
  sink->v[0] += sink->gpart->a_grav[0] * dt_kick_grav;
  sink->v[1] += sink->gpart->a_grav[1] * dt_kick_grav;
  sink->v[2] += sink->gpart->a_grav[2] * dt_kick_grav;

  /* Kick particles in momentum space (mesh forces) */
  sink->v[0] += sink->gpart->a_grav_mesh[0] * dt_kick_mesh_grav;
  sink->v[1] += sink->gpart->a_grav_mesh[1] * dt_kick_mesh_grav;
  sink->v[2] += sink->gpart->a_grav_mesh[2] * dt_kick_mesh_grav;

  /* Give the gpart friend the same velocity */
  sink->gpart->v_full[0] = sink->v[0];
  sink->gpart->v_full[1] = sink->v[1];
  sink->gpart->v_full[2] = sink->v[2];

  /* Kick extra variables */
  sink_kick_extra(sink, dt_kick_grav);
}

#endif /* SWIFT_KICK_H */
