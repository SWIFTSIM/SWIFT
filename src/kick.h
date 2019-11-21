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
#include "black_holes.h"
#include "const.h"
#include "debug.h"
#include "stars.h"
#include "timeline.h"
#include "boundary_movement.h"

/**
 * @brief Perform the 'kick' operation on a #gpart
 *
 * @param gp The #gpart to kick.
 * @param dt_kick_grav The kick time-step for gravity accelerations.
 * @param ti_start The starting (integer) time of the kick (for debugging
 * checks).
 * @param ti_end The ending (integer) time of the kick (for debugging checks).
 */
__attribute__((always_inline)) INLINE static void kick_gpart(
    struct gpart *restrict gp, double dt_kick_grav, integertime_t ti_start,
    integertime_t ti_end) {

#ifdef SWIFT_DEBUG_CHECKS
  if (gp->ti_kick != ti_start)
    error(
        "g-particle has not been kicked to the current time gp->ti_kick=%lld, "
        "ti_start=%lld, ti_end=%lld id=%lld",
        gp->ti_kick, ti_start, ti_end, gp->id_or_neg_offset);

  gp->ti_kick = ti_end;
#endif

  /* Kick particles in momentum space */
  gp->v_full[0] += gp->a_grav[0] * dt_kick_grav;
  gp->v_full[1] += gp->a_grav[1] * dt_kick_grav;
  gp->v_full[2] += gp->a_grav[2] * dt_kick_grav;

  /* Kick extra variables */
  gravity_kick_extra(gp, dt_kick_grav);
}

#ifdef WITH_ENGINEERING


/**
 * @brief Perform the 'kick' operation on a #part
 *
 * @param p The #part to kick.
 * @param xp The #xpart of the particle.
 * @param dt_kick_hydro The kick time-step for hydro accelerations.
 * @param dt_kick_grav The kick time-step for gravity accelerations.
 * @param dt_kick_therm The kick time-step for changes in thermal state.
 * @param dt_kick_corr The kick time-step for the gizmo-mfv gravity correction.
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme.
 * @param entropy_floor_props Properties of the entropy floor.
 * @param ti_start The starting (integer) time of the kick (for debugging
 * checks).
 * @param ti_end The ending (integer) time of the kick (for debugging checks).
 */
__attribute__((always_inline)) INLINE static void kick_part(
    struct part *restrict p, struct xpart *restrict xp, double dt_kick_hydro,
    double dt_kick_grav, double dt_kick_therm, double dt_kick_corr,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *entropy_floor_props,
    integertime_t ti_start, integertime_t ti_end) {
    if(p->since_euler >= 50) p->since_euler = 0;
    apply_boundary_motion(p, ti_start);
    if(!p->is_boundary){
      compute_particle_shifting(p, dt_kick_hydro);
    }

    if(p->since_euler == 0){
      double temp = p->rho;
      p->rho = p->rho + p->drho_dt*dt_kick_hydro;
      p->rho_t_minus1 = temp;
/*    if(p->id == 65744){
         printf("p->rho = %e, drho = %e, diff = %e\n", p->rho, p->drho_dt * dt_kick_hydro, p->rho - 1e-3);
      }*/
/*      if(!p->is_boundary){*/
        p->x[0] = p->x[0] + p->v[0]*dt_kick_hydro + 0.5*p->a_hydro[0]*dt_kick_hydro*dt_kick_hydro;
        p->x[1] = p->x[1] + p->v[1]*dt_kick_hydro + 0.5*p->a_hydro[1]*dt_kick_hydro*dt_kick_hydro;
        p->x[2] = p->x[2] + p->v[2]*dt_kick_hydro + 0.5*p->a_hydro[2]*dt_kick_hydro*dt_kick_hydro;
/*      }*/

      temp = p->v[0];
      p->v[0] = p->v[0] + p->a_hydro[0]*dt_kick_hydro;
      p->v_minus1[0] = temp;

      temp = p->v[1];
      p->v[1] = p->v[1] + p->a_hydro[1]*dt_kick_hydro;
      p->v_minus1[1] = temp;

      temp = p->v[2];
      p->v[2] = p->v[2] + p->a_hydro[2]*dt_kick_hydro;
      p->v_minus1[2] = temp;
#if defined(EOS_MULTIFLUID_TAIT)
      p->pressure = pressure_from_density(p->rho, p->rho_base);
/*  if(p->id == 462){
    printf("rho=%e, rho_base=%e, pressure=%e drho=%e\n", p->rho, p->rho_base, p->pressure, p->drho_dt*dt_kick_hydro);
  }*/
#else
      p->pressure = pressure_from_density(p->rho);
#endif


    }else{
      double temp = p->rho;
      p->rho = p->rho_t_minus1 + 2*dt_kick_hydro*p->drho_dt;
      p->rho_t_minus1 = temp;
//      if(!p->is_boundary){
        p->x[0] = p->x[0] + p->v[0]*dt_kick_hydro + 0.5*p->a_hydro[0]*dt_kick_hydro*dt_kick_hydro;
        p->x[1] = p->x[1] + p->v[1]*dt_kick_hydro + 0.5*p->a_hydro[1]*dt_kick_hydro*dt_kick_hydro;
        p->x[2] = p->x[2] + p->v[2]*dt_kick_hydro + 0.5*p->a_hydro[2]*dt_kick_hydro*dt_kick_hydro;
//      }
      temp = p->v[0];
      p->v[0] = p->v_minus1[0] + 2*p->a_hydro[0]*dt_kick_hydro;
      p->v_minus1[0] = temp;

      temp = p->v[1];
      p->v[1] = p->v_minus1[1] + 2*p->a_hydro[1]*dt_kick_hydro;
      p->v_minus1[1] = temp;

      temp = p->v[2];
      p->v[2] = p->v_minus1[2] + 2*p->a_hydro[2]*dt_kick_hydro;
      p->v_minus1[2] = temp;
#if defined(EOS_MULTIFLUID_TAIT)
      p->pressure = pressure_from_density(p->rho, p->rho_base);
  if(p->id == 462){
    printf("rho=%e, rho_base=%e, pressure=%e drho=%e\n", p->rho, p->rho_base, p->pressure, p->drho_dt*dt_kick_hydro);
  }
#else
      p->pressure = pressure_from_density(p->rho);
#endif

//      if(p->id == 65744){
//         printf("p->rho = %.11e, drho = %.11e, drho_dt = %.11e\n", p->rho, p->drho_dt * dt_kick_hydro, p->drho_dt);
//      }

    }
    p->since_euler += 1;
  /* Compute offsets since last cell construction */
  for (int k = 0; k < 3; k++) {
    const float dx = p->v[k]*dt_kick_hydro + 0.5*p->a_hydro[k]*dt_kick_hydro*dt_kick_hydro;
    if(dx != dx){
      error("NaN found %e %e %e %e %e %i", p->v[k], p->a_hydro[k], p->rho, p->rho_base, p->pressure, p->is_boundary);
    }
    xp->x_diff[k] -= dx;
    xp->x_diff_sort[k] -= dx;
  }


}
#else


/**
 * @brief Perform the 'kick' operation on a #part
 *
 * @param p The #part to kick.
 * @param xp The #xpart of the particle.
 * @param dt_kick_hydro The kick time-step for hydro accelerations.
 * @param dt_kick_grav The kick time-step for gravity accelerations.
 * @param dt_kick_therm The kick time-step for changes in thermal state.
 * @param dt_kick_corr The kick time-step for the gizmo-mfv gravity correction.
 * @param cosmo The cosmological model.
 * @param hydro_props The constants used in the scheme.
 * @param floor_props Properties of the entropy floor.
 * @param ti_start The starting (integer) time of the kick (for debugging
 * checks).
 * @param ti_end The ending (integer) time of the kick (for debugging checks).
 */
__attribute__((always_inline)) INLINE static void kick_part(
    struct part *restrict p, struct xpart *restrict xp, double dt_kick_hydro,
    double dt_kick_grav, double dt_kick_therm, double dt_kick_corr,
    const struct cosmology *cosmo, const struct hydro_props *hydro_props,
    const struct entropy_floor_properties *floor_props, integertime_t ti_start,
    integertime_t ti_end) {

#ifdef SWIFT_DEBUG_CHECKS
  if (p->ti_kick != ti_start)
    error(
        "particle has not been kicked to the current time p->ti_kick=%lld, "
        "ti_start=%lld, ti_end=%lld id=%lld time_bin=%d wakeup=%d",
        p->ti_kick, ti_start, ti_end, p->id, p->time_bin,
        p->limiter_data.wakeup);

  p->ti_kick = ti_end;
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

  /* Give the gpart friend the same velocity */
  if (p->gpart != NULL) {
    p->gpart->v_full[0] = xp->v_full[0];
    p->gpart->v_full[1] = xp->v_full[1];
    p->gpart->v_full[2] = xp->v_full[2];
  }

  /* Extra kick work */
  hydro_kick_extra(p, xp, dt_kick_therm, dt_kick_grav, dt_kick_hydro,
                   dt_kick_corr, cosmo, hydro_props, floor_props);
  if (p->gpart != NULL) gravity_kick_extra(p->gpart, dt_kick_grav);
#if !defined(EULER_ENG_SPH)
  /* Verify that the particle is not below the entropy floor */
  const float floor = entropy_floor(p, cosmo, entropy_floor_props);
  if (hydro_get_physical_entropy(p, xp, cosmo) < floor) {
    hydro_set_physical_entropy(p, xp, cosmo, floor);
    hydro_set_physical_internal_energy_dt(p, cosmo, 0.f);
  }
#endif
}
#endif
/**
 * @brief Perform the 'kick' operation on a #spart
 *
 * @param sp The #spart to kick.
 * @param dt_kick_grav The kick time-step for gravity accelerations.
 * @param ti_start The starting (integer) time of the kick (for debugging
 * checks).
 * @param ti_end The ending (integer) time of the kick (for debugging checks).
 */
__attribute__((always_inline)) INLINE static void kick_spart(
    struct spart *restrict sp, double dt_kick_grav, integertime_t ti_start,
    integertime_t ti_end) {

#ifdef SWIFT_DEBUG_CHECKS
  if (sp->ti_kick != ti_start)
    error(
        "s-particle has not been kicked to the current time sp->ti_kick=%lld, "
        "ti_start=%lld, ti_end=%lld id=%lld",
        sp->ti_kick, ti_start, ti_end, sp->id);

  sp->ti_kick = ti_end;
#endif

  /* Kick particles in momentum space */
  sp->v[0] += sp->gpart->a_grav[0] * dt_kick_grav;
  sp->v[1] += sp->gpart->a_grav[1] * dt_kick_grav;
  sp->v[2] += sp->gpart->a_grav[2] * dt_kick_grav;

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
 */
__attribute__((always_inline)) INLINE static void kick_bpart(
    struct bpart *restrict bp, double dt_kick_grav, integertime_t ti_start,
    integertime_t ti_end) {

#ifdef SWIFT_DEBUG_CHECKS
  if (bp->ti_kick != ti_start)
    error(
        "s-particle has not been kicked to the current time bp->ti_kick=%lld, "
        "ti_start=%lld, ti_end=%lld id=%lld",
        bp->ti_kick, ti_start, ti_end, bp->id);

  bp->ti_kick = ti_end;
#endif

  /* Kick particles in momentum space */
  bp->v[0] += bp->gpart->a_grav[0] * dt_kick_grav;
  bp->v[1] += bp->gpart->a_grav[1] * dt_kick_grav;
  bp->v[2] += bp->gpart->a_grav[2] * dt_kick_grav;

  /* Give the gpart friend the same velocity */
  bp->gpart->v_full[0] = bp->v[0];
  bp->gpart->v_full[1] = bp->v[1];
  bp->gpart->v_full[2] = bp->v[2];

  /* Kick extra variables */
  black_holes_kick_extra(bp, dt_kick_grav);
}

#endif /* SWIFT_KICK_H */
