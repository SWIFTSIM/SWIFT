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

    /* We want to make the particle finish it's time-step now. */

    /* Start by recovering the start and end point of the particle's time-step. */
    const integertime_t old_ti_beg =
            get_integer_time_begin(ti_current, p->time_bin);
    const integertime_t old_ti_end =
            get_integer_time_end(ti_current, p->time_bin);

    /* Old time-step length on the time-line */
    const integertime_t old_dti = old_ti_end - old_ti_beg;
    const integertime_t ti_end_half_old = old_ti_beg + old_dti / 2;

    /* The actual time-step size this particle will use */
    const integertime_t new_ti_beg = old_ti_beg;
    const integertime_t new_ti_end = ti_current;

#ifdef SWIFT_DEBUG_CHECKS
    const integertime_t new_dti = new_ti_end - new_ti_beg;

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

    /* Now we need to reverse the kick1...
     * Note the minus sign! (the dt are negative here) */
    dt_kick_hydro = -kick_get_hydro_kick_dt(old_ti_beg, ti_end_half_old,
                                            time_base, with_cosmology, cosmo);
    dt_kick_grav = -kick_get_grav_kick_dt(old_ti_beg, ti_end_half_old, time_base,
                                          with_cosmology, cosmo);
    dt_kick_therm = -kick_get_therm_kick_dt(old_ti_beg, ti_end_half_old,
                                            time_base, with_cosmology, cosmo);
    dt_kick_corr = -kick_get_corr_kick_dt(old_ti_beg, ti_end_half_old, time_base,
                                          with_cosmology, cosmo);

    /* Note that there is no need to change the mesh integration as we
     * can't go back more than one global step */
    kick_part(p, xp, dt_kick_hydro, dt_kick_grav, /*dt_kick_mesh_grav=*/0.,
              dt_kick_therm, dt_kick_corr, e->cosmology, e->hydro_properties,
              e->entropy_floor, ti_end_half_old, old_ti_beg,
            /*ti_start_mesh=*/-1, /*ti_end_mesh=*/-1);

    /* We can now produce a kick to the current point */
    dt_kick_hydro = kick_get_hydro_kick_dt(new_ti_beg, new_ti_end, time_base,
                                           with_cosmology, cosmo);
    dt_kick_grav = kick_get_grav_kick_dt(new_ti_beg, new_ti_end, time_base,
                                         with_cosmology, cosmo);
    dt_kick_therm = kick_get_therm_kick_dt(new_ti_beg, new_ti_end, time_base,
                                           with_cosmology, cosmo);
    dt_kick_corr = kick_get_corr_kick_dt(new_ti_beg, new_ti_end, time_base,
                                         with_cosmology, cosmo);

    /* Note that there is no need to change the mesh integration as we
     * can't go back more than one global step */
    kick_part(p, xp, dt_kick_hydro, dt_kick_grav, /*dt_kick_mesh_grav=*/0.,
              dt_kick_therm, dt_kick_corr, e->cosmology, e->hydro_properties,
              e->entropy_floor, new_ti_beg, new_ti_end,
            /*ti_start_mesh=*/-1, /*ti_end_mesh=*/-1);

    /* The particle is now ready to compute its new time-step size and for the
     * next kick */
    p->time_bin = -min_active_bin;
    p->limiter_data.wakeup = time_bin_not_awake;
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
INLINE static void timestep_process_sync_dmpart(struct dmpart *p, const struct engine *e,
                                                const struct cosmology *cosmo) {

    const timebin_t min_active_bin = e->min_active_bin;
    const timebin_t max_active_bin = e->max_active_bin;

    /* This particle is already active. Nothing to do here... */
    if (p->time_bin <= max_active_bin) {
        return;
    }

    /* Did this particle had a SIDM kick? if so, resolve and skip kick2 */
    sidm_kick_to_dmpart(p);

#ifdef SWIFT_DEBUG_CHECKS
    /* With SIDM kick particle has been kicked to current time */
  const integertime_t ti_current = e->ti_current;
  const integertime_t ti_end = ti_current;
  p->ti_kick = ti_end;
#endif


    /* The particle is now ready to compute its new time-step size and for the
     * next kick */
    p->time_bin = -min_active_bin;
    p->limiter_data.wakeup = time_bin_not_awake;

}


#endif /* SWIFT_TIMESTEP_SYNC_H */
