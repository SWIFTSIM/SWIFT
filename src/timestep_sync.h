/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2018 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#ifndef SWIFT_TIMESTEP_SYNC_H
#define SWIFT_TIMESTEP_SYNC_H

/* Config parameters. */
#include "../config.h"

/* Local includes */
#include "engine.h"
#include "kick.h"
#include "timestep.h"
#include "drift.h"

/**
 * @brief Processes a particle that has been flagged for synchronization on the
 * time-line.
 *
 * We revert the particle's kick and apply a new one that ends at the current
 * time. The particle is then ready to compute a new time-step and proceed with
 * a regular kick1.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param e The #engine.
 * @param cosmo The cosmology model.
 */
INLINE static void timestep_process_sync_dmpart(struct dmpart *p, const struct engine *e,
                                              const struct cosmology *cosmo) {

  const int with_cosmology = (e->policy & engine_policy_cosmology);
  const integertime_t ti_current = e->ti_current;
  const timebin_t min_active_bin = e->min_active_bin;
  const timebin_t max_active_bin = e->max_active_bin;
  const double time_base = e->time_base;
    
  /* This particle is already active. Nothing to do here... */
  if (p->time_bin <= max_active_bin) {
      return;
  }

  /* We want to make the particle finish it's time-step now. */

  /* Start by recovering the start and end point of the particle's time-step. */
  const integertime_t old_ti_beg =
      get_integer_time_begin(ti_current, p->time_bin);
  const integertime_t old_ti_end =
      get_integer_time_end(ti_current, p->time_bin);

  /* Old time-step length on the time-line */
  const integertime_t old_dti = old_ti_end - old_ti_beg;

  /* The actual time-step size this particle will use */
  const integertime_t new_ti_beg = old_ti_beg;
  const integertime_t new_ti_end = ti_current;
  const integertime_t new_dti = new_ti_end - new_ti_beg;

#ifdef SWIFT_DEBUG_CHECKS
  /* Some basic safety checks */
  if (old_ti_beg >= ti_current)
    error(
        "Incorrect value for old time-step beginning ti_current=%lld, "
        "old_ti_beg=%lld",
        ti_current, old_ti_beg);

  if (old_ti_end <= ti_current)
    error(
        "Incorrect value for old time-step end ti_current=%lld, "
        "old_ti_end=%lld",
        ti_current, old_ti_end);

  if (new_ti_end > old_ti_end) error("New end of time-step after the old one");

  if (new_dti > old_dti) error("New time-step larger than old one");
#endif

  double dt_kick_grav = 0.;

  /* Now we need to reverse the kick1... (the dt are negative here) */
  if (with_cosmology) {
    dt_kick_grav = -cosmology_get_grav_kick_factor(cosmo, old_ti_beg,
                                                   old_ti_beg + old_dti / 2);
  } else {
    dt_kick_grav = -(old_dti / 2) * time_base;
  }

  kick_dmpart(p, dt_kick_grav, old_ti_beg + old_dti / 2, old_ti_beg);

  /* We can now produce a kick to the current point */
  if (with_cosmology) {
    dt_kick_grav =
        cosmology_get_grav_kick_factor(cosmo, old_ti_beg, new_ti_beg + new_dti);
      
  } else {
    dt_kick_grav = (new_dti) * time_base;
  }

  kick_dmpart(p, dt_kick_grav, new_ti_beg, new_ti_beg + new_dti);
  
  double dt_drift = 0.;
    
  if (with_cosmology) {
      dt_drift = cosmology_get_drift_factor(cosmo, old_ti_beg, new_ti_beg + new_dti);
  } else {
      dt_drift = (new_dti) * e->time_base;
  }
    
  /* Drift dm part to current time */
  /*drift_dmpart(p, dt_drift, old_ti_beg, new_ti_end);*/
    
  /* Did this particle had a SIDM kick? if so, resolve */
  do_sidm_kick_to_dmpart(p, dt_drift);

  /* The particle is now ready to compute its new time-step size and for the
   * next kick */
  p->time_bin = -min_active_bin;
  p->limiter_data.wakeup = time_bin_not_awake;
  p->limiter_data.to_be_synchronized = 0;

}


/**
 * @brief Processes a particle that has been flagged for synchronization on the
 * time-line.
 *
 * We revert the particle's kick and apply a new one that ends at the current
 * time. The particle is then ready to compute a new time-step and proceed with
 * a regular kick1.
 *
 * @param p The #part.
 * @param xp The #xpart.
 * @param e The #engine.
 * @param cosmo The cosmology model.
 */
INLINE static void timestep_process_sync_part(struct part *p, struct xpart *xp,
                                              const struct engine *e,
                                              const struct cosmology *cosmo) {
    
    const int with_cosmology = (e->policy & engine_policy_cosmology);
    const integertime_t ti_current = e->ti_current;
    const timebin_t max_active_bin = e->max_active_bin;
    const timebin_t min_active_bin = e->min_active_bin;
    const double time_base = e->time_base;
    
    p->limiter_data.to_be_synchronized = 0;
    
    /* This particle is already active. Nothing to do here... */
    if (p->time_bin <= max_active_bin) {
        return;
    }
    
    // message(" Synchronizing particle! %lld old bin=%d", p->id, p->time_bin);
    
    /* We want to make the particle finish it's time-step now. */
    
    /* Start by recovering the start and end point of the particle's time-step. */
    const integertime_t old_ti_beg =
    get_integer_time_begin(ti_current, p->time_bin);
    const integertime_t old_ti_end =
    get_integer_time_end(ti_current, p->time_bin);
    
    /* Old time-step length on the time-line */
    const integertime_t old_dti = old_ti_end - old_ti_beg;
    
    /* The actual time-step size this particle will use */
    const integertime_t new_ti_beg = old_ti_beg;
    const integertime_t new_ti_end = ti_current;
    const integertime_t new_dti = new_ti_end - new_ti_beg;
    
#ifdef SWIFT_DEBUG_CHECKS
    /* Some basic safety checks */
    if (old_ti_beg >= ti_current)
        error(
              "Incorrect value for old time-step beginning ti_current=%lld, "
              "old_ti_beg=%lld",
              ti_current, old_ti_beg);
    
    if (old_ti_end <= ti_current)
        error(
              "Incorrect value for old time-step end ti_current=%lld, "
              "old_ti_end=%lld",
              ti_current, old_ti_end);
    
    if (new_ti_end > old_ti_end) error("New end of time-step after the old one");
    
    if (new_dti > old_dti) error("New time-step larger than old one");
#endif
    
    double dt_kick_grav = 0., dt_kick_hydro = 0., dt_kick_therm = 0.,
    dt_kick_corr = 0.;
    
    /* Now we need to reverse the kick1... (the dt are negative here) */
    if (with_cosmology) {
        dt_kick_hydro = -cosmology_get_hydro_kick_factor(cosmo, old_ti_beg,
                                                         old_ti_beg + old_dti / 2);
        dt_kick_grav = -cosmology_get_grav_kick_factor(cosmo, old_ti_beg,
                                                       old_ti_beg + old_dti / 2);
        dt_kick_therm = -cosmology_get_therm_kick_factor(cosmo, old_ti_beg,
                                                         old_ti_beg + old_dti / 2);
        dt_kick_corr = -cosmology_get_corr_kick_factor(cosmo, old_ti_beg,
                                                       old_ti_beg + old_dti / 2);
    } else {
        dt_kick_hydro = -(old_dti / 2) * time_base;
        dt_kick_grav = -(old_dti / 2) * time_base;
        dt_kick_therm = -(old_dti / 2) * time_base;
        dt_kick_corr = -(old_dti / 2) * time_base;
    }
    
    kick_part(p, xp, dt_kick_hydro, dt_kick_grav, dt_kick_therm, dt_kick_corr,
              e->cosmology, e->hydro_properties, e->entropy_floor,
              old_ti_beg + old_dti / 2, old_ti_beg);
    
    /* We can now produce a kick to the current point */
    if (with_cosmology) {
        dt_kick_hydro = cosmology_get_hydro_kick_factor(cosmo, new_ti_beg,
                                                        new_ti_beg + new_dti);
        dt_kick_grav =
        cosmology_get_grav_kick_factor(cosmo, old_ti_beg, new_ti_beg + new_dti);
        dt_kick_therm = cosmology_get_therm_kick_factor(cosmo, old_ti_beg,
                                                        new_ti_beg + new_dti);
        dt_kick_corr =
        cosmology_get_corr_kick_factor(cosmo, old_ti_beg, new_ti_beg + new_dti);
    } else {
        dt_kick_hydro = (new_dti)*time_base;
        dt_kick_grav = (new_dti)*time_base;
        dt_kick_therm = (new_dti)*time_base;
        dt_kick_corr = (new_dti)*time_base;
    }
    
    kick_part(p, xp, dt_kick_hydro, dt_kick_grav, dt_kick_therm, dt_kick_corr,
              e->cosmology, e->hydro_properties, e->entropy_floor, new_ti_beg,
              new_ti_beg + new_dti);
    
    /* The particle is now ready to compute its new time-step size and for the
     * next kick */
    p->time_bin = -min_active_bin;
    p->limiter_data.wakeup = time_bin_not_awake;
}

#endif /* SWIFT_TIMESTEP_SYNC_H */
