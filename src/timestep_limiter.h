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
#ifndef SWIFT_TIMESTEP_LIMITER_H
#define SWIFT_TIMESTEP_LIMITER_H

/* Config parameters. */
#include "../config.h"

/* Local headers. */
#include "kick.h"

/**
 * @brief Prepare the limiter-quantities for a sidm calculation.
 *
 * @param p The #dmpart.
 */
__attribute__((always_inline)) INLINE static void
timestep_limiter_prepare_sidm(struct dmpart *restrict p) {
    
    p->limiter_data.min_ngb_time_bin = num_time_bins + 1;
}

/**
 * @brief Prepare the limiter-quantities for a hydro force calculation.
 *
 * @param p The #part.
 * @param xp The #xpart.
 */
__attribute__((always_inline)) INLINE static void
timestep_limiter_prepare_force(struct part *restrict p,
                               struct xpart *restrict xp) {

  p->limiter_data.min_ngb_time_bin = num_time_bins + 1;
}

/**
 * @brief Terminate the limiter-quantities after a hydro force calculation.
 *
 * @param p The #part.
 */
__attribute__((always_inline)) INLINE static void timestep_limiter_end_force(
    struct part *restrict p) {
#ifdef SWIFT_DEBUG_CHECKS
  if (p->limiter_data.min_ngb_time_bin == 0)
    error("Minimal time-bin of neighbours is 0");
#endif
}

/**
 * @brief Wakes up a particle by rewinding it's kick1 back in time and applying
 * a new one such that the particle becomes active again in the next time-step.
 *
 * @param p The #part to update.
 * @param xp Its #xpart companion.
 * @param e The #engine (to extract time-line information).
 *
 * @return The updated integer end-of-step of the particle.
 */
__attribute__((always_inline)) INLINE static integertime_t timestep_limit_part(
    struct part *restrict p, struct xpart *restrict xp,
    const struct engine *e) {

  const struct cosmology *cosmo = e->cosmology;
  const int with_cosmology = e->policy & engine_policy_cosmology;
  const double time_base = e->time_base;
  const integertime_t ti_current = e->ti_current;

  if (part_is_active(p, e)) {

    /* First case, the particle was active so we only need to update the length
       of its next time-step */

    /* New time-bin of this particle */
    p->time_bin = -p->limiter_data.wakeup + 2;

    /* Mark the particle as being rady to be time integrated */
    p->limiter_data.wakeup = time_bin_not_awake;

    /* Return the new end-of-step for this particle */
    return ti_current + get_integer_timestep(p->time_bin);

  } else {

    /* Second case, the particle was inactive so we need to interrupt its
       time-step, undo the "kick" operator and assign a new time-step size */

    /* The timebins to play with */
    const timebin_t old_bin = p->time_bin;
    const timebin_t new_bin = -p->limiter_data.wakeup + 2;

    /* Current start and end time of this particle */
    const integertime_t ti_beg_old =
        get_integer_time_begin(ti_current, old_bin);
    const integertime_t ti_end_old = get_integer_time_end(ti_current, old_bin);

    /* Length of the old and new time-step */
    const integertime_t dti_old = ti_end_old - ti_beg_old;
    const integertime_t dti_new = get_integer_timestep(new_bin);

    /* Let's now search for the starting point of the new step */
    int k = 0;
    while (ti_beg_old + k * dti_new <= ti_current) k++;

    const integertime_t ti_beg_new = ti_beg_old + (k - 1) * dti_new;

    /* End of the old half step */
    const integertime_t ti_end_half_old = ti_beg_old + dti_old / 2;

    /* End of the new half step */
    const integertime_t ti_end_half_new = ti_beg_new + dti_new / 2;

#ifdef SWIFT_DEBUG_CHECKS
    /* Some basic safety checks */
    if (ti_beg_old >= e->ti_current)
      error(
          "Incorrect value for old time-step beginning ti_current=%lld, "
          "ti_beg_old=%lld",
          e->ti_current, ti_beg_old);

    if (ti_end_old <= e->ti_current)
      error(
          "Incorrect value for old time-step end ti_current=%lld, "
          "ti_end_old=%lld",
          e->ti_current, ti_end_old);

    if (ti_beg_new < ti_beg_old)
      error("New beg of time-step before the old one");

    if (dti_new > dti_old) error("New time-step larger than old one");
#endif

    double dt_kick_grav = 0., dt_kick_hydro = 0., dt_kick_therm = 0.,
           dt_kick_corr = 0.;

    /* Now we need to reverse the kick1...
     * Note the minus sign! (the dt are negative here) */
    dt_kick_hydro = -kick_get_hydro_kick_dt(ti_beg_old, ti_end_half_old,
                                            time_base, with_cosmology, cosmo);
    dt_kick_grav = -kick_get_grav_kick_dt(ti_beg_old, ti_end_half_old,
                                          time_base, with_cosmology, cosmo);
    dt_kick_therm = -kick_get_therm_kick_dt(ti_beg_old, ti_end_half_old,
                                            time_base, with_cosmology, cosmo);
    dt_kick_corr = -kick_get_corr_kick_dt(ti_beg_old, ti_end_half_old,
                                          time_base, with_cosmology, cosmo);

    /* Note that there is no need to change the mesh integration as we
     * can't go back more than one global step */
    kick_part(p, xp, dt_kick_hydro, dt_kick_grav, /*dt_kick_mesh_grav=*/0.,
              dt_kick_therm, dt_kick_corr, e->cosmology, e->hydro_properties,
              e->entropy_floor, ti_end_half_old, ti_beg_old,
              /*ti_start_mesh=*/-1, /*ti_end_mesh=*/-1);

    /* ...and apply the new one (dt is positiive).
     * This brings us to the current time. */
    dt_kick_hydro = kick_get_hydro_kick_dt(ti_beg_old, ti_beg_new, time_base,
                                           with_cosmology, cosmo);
    dt_kick_grav = kick_get_grav_kick_dt(ti_beg_old, ti_beg_new, time_base,
                                         with_cosmology, cosmo);
    dt_kick_therm = kick_get_therm_kick_dt(ti_beg_old, ti_beg_new, time_base,
                                           with_cosmology, cosmo);
    dt_kick_corr = kick_get_corr_kick_dt(ti_beg_old, ti_beg_new, time_base,
                                         with_cosmology, cosmo);

    /* Note that there is no need to change the mesh integration as we
     * can't go back more than one global step */
    kick_part(p, xp, dt_kick_hydro, dt_kick_grav, /*dt_kick_mesh_grav=*/0.,
              dt_kick_therm, dt_kick_corr, e->cosmology, e->hydro_properties,
              e->entropy_floor, ti_beg_old, ti_beg_new, /*ti_start_mesh=*/-1,
              /*ti_end_mesh=*/-1);

    /* The particle has now been kicked to the current time */

    /* New time-bin of this particle */
    p->time_bin = new_bin;

    /* Mark the particle as being ready to be time integrated */
    p->limiter_data.wakeup = time_bin_not_awake;

    /* Do we need to apply the mising kick1 or is this bin active and
       it will be done in the kick task? */
    if (new_bin > e->max_active_bin) {

      /* Apply the missing kick1 */

      dt_kick_hydro = kick_get_hydro_kick_dt(ti_beg_new, ti_end_half_new,
                                             time_base, with_cosmology, cosmo);
      dt_kick_grav = kick_get_grav_kick_dt(ti_beg_new, ti_end_half_new,
                                           time_base, with_cosmology, cosmo);
      dt_kick_therm = kick_get_therm_kick_dt(ti_beg_new, ti_end_half_new,
                                             time_base, with_cosmology, cosmo);
      dt_kick_corr = kick_get_corr_kick_dt(ti_beg_new, ti_end_half_new,
                                           time_base, with_cosmology, cosmo);

      /* Note that there is no need to change the mesh integration as we
       * can't go back more than one global step */
      kick_part(p, xp, dt_kick_hydro, dt_kick_grav, /*dt_kick_mesh_grav=*/0.,
                dt_kick_therm, dt_kick_corr, e->cosmology, e->hydro_properties,
                e->entropy_floor, ti_beg_new, ti_end_half_new,
                /*ti_start_mesh=*/-1, /*ti_end_mesh=*/-1);

      /* Return the new end-of-step for this particle */
      return ti_beg_new + dti_new;

    } else {

      /* No kick to do here */

      /* Return the new end-of-step for this particle */
      return ti_current + get_integer_timestep(p->time_bin);
    }
  }
}

/**
 * @brief Wakes up a particle by rewinding it's kick1 back in time and applying
 * a new one such that the particle becomes active again in the next time-step.
 *
 * @param p The #part to update.
 * @param xp Its #xpart companion.
 * @param e The #engine (to extract time-line information).
 *
 * @return The updated integer end-of-step of the particle.
 */
__attribute__((always_inline)) INLINE static integertime_t timestep_limit_dmpart(
    struct dmpart *restrict p, const struct engine *e) {
    
    const struct cosmology *cosmo = e->cosmology;
    const int with_cosmology = e->policy & engine_policy_cosmology;
    const double time_base = e->time_base;
    const integertime_t ti_current = e->ti_current;
    
    if (dmpart_is_active(p, e)) {
        
        /* First case, the particle was active so we only need to update the length
         of its next time-step */
        
        /* New time-bin of this particle */
        p->time_bin = -p->limiter_data.wakeup + 2;
        
        /* Mark the particle as being rady to be time integrated */
        p->limiter_data.wakeup = time_bin_not_awake;
        
        /* Return the new end-of-step for this particle */
        return ti_current + get_integer_timestep(p->time_bin);
        
    } else {
        
        /* Second case, the particle was inactive so we need to interrupt its
         time-step, undo the "kick" operator and assign a new time-step size */
        
        /* The timebins to play with */
        const timebin_t old_bin = p->time_bin;
        const timebin_t new_bin = -p->limiter_data.wakeup + 2;
        
        /* Current start and end time of this particle */
        const integertime_t ti_beg_old = get_integer_time_begin(ti_current, old_bin);
        const integertime_t ti_end_old = get_integer_time_end(ti_current, old_bin);
        
        /* Length of the old and new time-step */
        const integertime_t dti_old = ti_end_old - ti_beg_old;
        const integertime_t dti_new = get_integer_timestep(new_bin);
        
        /* Let's now search for the starting point of the new step */
        int k = 0;
        while (ti_beg_old + k * dti_new <= ti_current) k++;
        
        const integertime_t ti_beg_new = ti_beg_old + (k - 1) * dti_new;
        
        double dt_kick_grav = 0.;
        
        /* Now we need to reverse the kick1... (the dt are negative here) */
        if (with_cosmology) {
            dt_kick_grav = -cosmology_get_grav_kick_factor(cosmo, ti_beg_old, ti_beg_old + dti_old / 2);
        } else {
            dt_kick_grav = -(dti_old / 2) * time_base;
        }
        
        kick_dmpart(p, dt_kick_grav, ti_beg_old + dti_old / 2, ti_beg_old);
        
        /* ...and apply the new one (dt is positive).
         * This brings us to the current time. */
        if (with_cosmology) {
            dt_kick_grav = cosmology_get_grav_kick_factor(cosmo, ti_beg_new, ti_beg_new + dti_new / 2);
        } else {
            dt_kick_grav = (ti_beg_new - ti_beg_old) * time_base;
        }
        
        kick_dmpart(p, dt_kick_grav, ti_beg_old, ti_beg_new);
        
        /* The particle has now been kicked to the current time */
        
        /* New time-bin of this particle */
        p->time_bin = new_bin;
        
        /* Mark the particle as being ready to be time integrated */
        p->limiter_data.wakeup = time_bin_not_awake;
        
        /* Do we need to apply the mising kick1 or is this bin active and
         it will be done in the kick task? */
        if (new_bin > e->max_active_bin) {
            
            /* Apply the missing kick1 */
            
            if (with_cosmology) {
                dt_kick_grav = cosmology_get_grav_kick_factor(cosmo, ti_beg_new, ti_beg_new + dti_new / 2);
            } else {
                dt_kick_grav = (dti_new / 2) * time_base;
            }
            
            kick_dmpart(p, dt_kick_grav, ti_beg_new, ti_beg_new + dti_new / 2);
            
            /* Return the new end-of-step for this particle */
            return ti_beg_new + dti_new;
            
        } else {
            
            /* No kick to do here */
            
            /* Return the new end-of-step for this particle */
            return ti_current + get_integer_timestep(p->time_bin);
        }
    }
}

#endif /* SWIFT_TIMESTEP_LIMITER_H */
