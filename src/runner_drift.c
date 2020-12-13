/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2012 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *                    Matthieu Schaller (matthieu.schaller@durham.ac.uk)
 *               2015 Peter W. Draper (p.w.draper@durham.ac.uk)
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

/* Config parameters. */
#include "../config.h"

/* This object's header. */
#include "runner.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "timers.h"
#include "dark_matter.h"
#include "drift.h"


/**
 * @brief Drift all part in a cell.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_drift_part(struct runner *r, struct cell *c, int timer) {

  TIMER_TIC;

  cell_drift_part(c, r->e, 0);

  if (timer) TIMER_TOC(timer_drift_part);
}

/**
 * @brief Drift all gpart in a cell.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_drift_gpart(struct runner *r, struct cell *c, int timer) {

  TIMER_TIC;

  cell_drift_gpart(c, r->e, 0);

  if (timer) TIMER_TOC(timer_drift_gpart);
}

/**
 * @brief Drift all spart in a cell.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_drift_spart(struct runner *r, struct cell *c, int timer) {

  TIMER_TIC;

  cell_drift_spart(c, r->e, 0);

  if (timer) TIMER_TOC(timer_drift_spart);
}

/**
 * @brief Drift all bpart in a cell.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_drift_bpart(struct runner *r, struct cell *c, int timer) {

  TIMER_TIC;

  cell_drift_bpart(c, r->e, 0);

  if (timer) TIMER_TOC(timer_drift_bpart);
}

/**
 * @brief Drift all dmpart in a cell.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_drift_dmpart(struct runner *r, struct cell *c, int timer) {
    
    /*TIMER_TIC;
    
    cell_drift_dmpart(c, r->e, 0);
    
    if (timer) TIMER_TOC(timer_drift_dmpart);*/
    
    TIMER_TIC;
    
    const struct engine *e = r->e;
    const int periodic = e->s->periodic;
    const double dim[3] = {e->s->dim[0], e->s->dim[1], e->s->dim[2]};
    const int with_cosmology = (e->policy & engine_policy_cosmology);
    const float dark_matter_h_max = e->sidm_properties->h_max;
    const float dark_matter_h_min = e->sidm_properties->h_min;
    const integertime_t ti_old_dmpart = c->dark_matter.ti_old_part;
    const integertime_t ti_current = e->ti_current;
    struct dmpart *const dmparts = c->dark_matter.parts;
    
    float dx_max = 0.f, dx2_max = 0.f;
    float cell_h_max = 0.f;
    
    /* Early abort? */
    if (c->dark_matter.count == 0) {
        
        /* Update the time of the last drift */
        /*c->dark_matter.ti_old_part = ti_current;*/
        
        return;
    }
    
    /* Ok, we have some particles somewhere in the hierarchy to drift */
    
    /* Are we not in a leaf ? */
    if (c->split) {
        
        /* Loop over the progeny and collect their data. */
        for (int k = 0; k < 8; k++) {
            if (c->progeny[k] != NULL) {
                struct cell *cp = c->progeny[k];
                
                /* Recurse */
                runner_do_drift_dmpart(r, cp, 1);
                
                /* Update */
                dx_max = max(dx_max, cp->dark_matter.dx_max_part);
                cell_h_max = max(cell_h_max, cp->dark_matter.h_max);
            }
        }
        
        /* Store the values */
        c->dark_matter.h_max = cell_h_max;
        c->dark_matter.dx_max_part = dx_max;
        
        /* Update the time of the last drift */
        /*c->dark_matter.ti_old_part = ti_current;*/

        
    } else if (!c->split && ti_current > ti_old_dmpart) {
        
        /* Drift from the last time the cell was drifted to the current time */
        double dt_drift;
        if (with_cosmology) {
            dt_drift =
            cosmology_get_drift_factor(e->cosmology, ti_old_dmpart, ti_current);
        } else {
            dt_drift = (ti_current - ti_old_dmpart) * e->time_base;
        }

        dt_drift /= 2.f;
        
        /* Loop over all the star particles in the cell */
        const size_t nr_dmparts = c->dark_matter.count;
        for (size_t k = 0; k < nr_dmparts; k++) {
            /* Get a handle on the spart. */
            struct dmpart *const dmp = &dmparts[k];
            
            /* Ignore inhibited particles */
            if (dmpart_is_inhibited(dmp, e)) continue;
            
            /* Drift... */
            drift_dmpart(dmp, dt_drift, ti_old_dmpart, ti_current);
            
#ifdef SWIFT_DEBUG_CHECKS
            /* Make sure the particle does not drift by more than a box length. */
            if (fabs(dmp->v_full[0] * dt_drift) > e->s->dim[0] ||
                fabs(dmp->v_full[1] * dt_drift) > e->s->dim[1] ||
                fabs(dmp->v_full[2] * dt_drift) > e->s->dim[2]) {
                error("DM Particle drifts by more than a box length!");
            }
#endif
            
            /* In non-periodic BC runs, remove particles that crossed the border */
            if (!periodic) {
                
                /* Did the particle leave the box?  */
                if ((dmp->x[0] > dim[0]) || (dmp->x[0] < 0.) ||  // x
                    (dmp->x[1] > dim[1]) || (dmp->x[1] < 0.) ||  // y
                    (dmp->x[2] > dim[2]) || (dmp->x[2] < 0.)) {  // z
                    
                    lock_lock(&e->s->lock);
                    
                    if (lock_unlock(&e->s->lock) != 0)
                    error("Failed to unlock the space!");
                    
                    continue;
                }
            }
            
            /* Limit h to within the allowed range */
            dmp->h = min(dmp->h, dark_matter_h_max);
            dmp->h = max(dmp->h, dark_matter_h_min);
            
            /* Compute (square of) motion since last cell construction */
            const float dx2 = dmp->x_diff[0] * dmp->x_diff[0] +
            dmp->x_diff[1] * dmp->x_diff[1] +
            dmp->x_diff[2] * dmp->x_diff[2];
            dx2_max = max(dx2_max, dx2);
            
            /* Maximal smoothing length */
            cell_h_max = max(cell_h_max, dmp->h);
            
            /* Get ready for a density calculation */
            if (dmpart_is_active(dmp, e)) dark_matter_init_dmpart(dmp);
            
            /* All dmparts get ready for a SIDM calculation */
            sidm_init_dmpart(dmp);
        }
        
        /* Now, get the maximal particle motion from its square */
        dx_max = sqrtf(dx2_max);
        
        /* Store the values */
        c->dark_matter.h_max = cell_h_max;
        c->dark_matter.dx_max_part = dx_max;
        
        /* Update the time of the last drift */
        /*c->dark_matter.ti_old_part = ti_current;*/

        
    }
    
    if (timer) TIMER_TOC(timer_drift_dmpart);

}


/**
 * @brief Drift all dmpart in a cell.
 *
 * @param r The runner thread.
 * @param c The cell.
 * @param timer Are we timing this ?
 */
void runner_do_drift2_dmpart(struct runner *r, struct cell *c) {
    
    const struct engine *e = r->e;
    const int periodic = e->s->periodic;
    const double dim[3] = {e->s->dim[0], e->s->dim[1], e->s->dim[2]};
    const int with_cosmology = (e->policy & engine_policy_cosmology);
    const float dark_matter_h_max = e->sidm_properties->h_max;
    const float dark_matter_h_min = e->sidm_properties->h_min;
    const integertime_t ti_old_dmpart = c->dark_matter.ti_old_part;
    const integertime_t ti_current = e->ti_current;
    struct dmpart *const dmparts = c->dark_matter.parts;
    
    float dx_max = 0.f, dx2_max = 0.f;
    float cell_h_max = 0.f;
    
    /* Early abort? */
    if (c->dark_matter.count == 0) {
        
        /* Update the time of the last drift */
        c->dark_matter.ti_old_part = ti_current;

        return;
    }
    
    /* Ok, we have some particles somewhere in the hierarchy to drift */
    
    /* Are we not in a leaf ? */
    if (c->split) {
        
        /* Loop over the progeny and collect their data. */
        for (int k = 0; k < 8; k++) {
            if (c->progeny[k] != NULL) {
                struct cell *cp = c->progeny[k];
                
                /* Recurse */
                runner_do_drift2_dmpart(r, cp);
                
                /* Update */
                dx_max = max(dx_max, cp->dark_matter.dx_max_part);
                cell_h_max = max(cell_h_max, cp->dark_matter.h_max);
            }
        }
        
        /* Store the values */
        c->dark_matter.h_max = cell_h_max;
        c->dark_matter.dx_max_part = dx_max;
        
        /* Update the time of the last drift */
        c->dark_matter.ti_old_part = ti_current;

        
    } else if (!c->split) {

        /* Drift from the last time the cell was drifted to the current time */
        double dt_drift;
        if (with_cosmology) {
            dt_drift =
            cosmology_get_drift_factor(e->cosmology, ti_old_dmpart, ti_current);
        } else {
            dt_drift = (ti_current - ti_old_dmpart) * e->time_base;
        }

        dt_drift /= 2.f;
        
        /* Loop over all the star particles in the cell */
        const size_t nr_dmparts = c->dark_matter.count;
        for (size_t k = 0; k < nr_dmparts; k++) {
            /* Get a handle on the spart. */
            struct dmpart *const dmp = &dmparts[k];
            
            /* Ignore inhibited particles */
            if (dmpart_is_inhibited(dmp, e)) continue;
            
            /* Drift... */
            drift_dmpart(dmp, dt_drift, ti_old_dmpart, ti_current);

#ifdef SWIFT_DEBUG_CHECKS
            /* Make sure the particle does not drift by more than a box length. */
            if (fabs(dmp->v_full[0] * dt_drift) > e->s->dim[0] ||
                fabs(dmp->v_full[1] * dt_drift) > e->s->dim[1] ||
                fabs(dmp->v_full[2] * dt_drift) > e->s->dim[2]) {
                error("DM Particle drifts by more than a box length!");
            }
#endif
            
            /* In non-periodic BC runs, remove particles that crossed the border */
            if (!periodic) {
                
                /* Did the particle leave the box?  */
                if ((dmp->x[0] > dim[0]) || (dmp->x[0] < 0.) ||  // x
                    (dmp->x[1] > dim[1]) || (dmp->x[1] < 0.) ||  // y
                    (dmp->x[2] > dim[2]) || (dmp->x[2] < 0.)) {  // z
                    
                    lock_lock(&e->s->lock);
                    
                    if (lock_unlock(&e->s->lock) != 0)
                    error("Failed to unlock the space!");
                    
                    continue;
                }
            }
            
            /* Limit h to within the allowed range */
            dmp->h = min(dmp->h, dark_matter_h_max);
            dmp->h = max(dmp->h, dark_matter_h_min);
            
            /* Compute (square of) motion since last cell construction */
            const float dx2 = dmp->x_diff[0] * dmp->x_diff[0] +
            dmp->x_diff[1] * dmp->x_diff[1] +
            dmp->x_diff[2] * dmp->x_diff[2];
            dx2_max = max(dx2_max, dx2);
            
            /* Maximal smoothing length */
            cell_h_max = max(cell_h_max, dmp->h);
        }
        
        /* Now, get the maximal particle motion from its square */
        dx_max = sqrtf(dx2_max);
        
        /* Store the values */
        c->dark_matter.h_max = cell_h_max;
        c->dark_matter.dx_max_part = dx_max;
        
        /* Update the time of the last drift */
        c->dark_matter.ti_old_part = ti_current;

    }
    
}



