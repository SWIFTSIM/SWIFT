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
#include "dark_matter.h"
#include "dark_matter_iact.h"
#include "cell.h"
#include "engine.h"
#include "timers.h"

/**
 * @brief
 * @param r The thread #runner.
 * @param c The #cell.
 * @param timer Are we timing this?
 */
void runner_do_self_dark_matter_SIDM(struct runner *r, struct cell *c, int timer) {
    
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    const int with_cosmology = e->policy & engine_policy_cosmology;
    const struct unit_system *us = e->internal_units;
    const struct sidm_props *sidm_props = e->sidm_properties;
    
    TIMER_TIC;
    
    /* Anything to do here? */
    if (!cell_is_active_dark_matter(c, e)) return;
    
    /* Cosmological terms */
    const float a = cosmo->a;
    const float H = cosmo->H;
    
    const int count = c->dark_matter.count;
    struct dmpart *restrict dmparts = c->dark_matter.parts;
    
    /* Loop over the parts in ci. */
    for (int pid = 0; pid < count; pid++) {
        
        /* Get a hold of the ith part in ci. */
        struct dmpart *restrict pi = &dmparts[pid];
        
        /* Skip inhibited particles. */
        if (dmpart_is_inhibited(pi, e)) continue;
        
        const int pi_active = dmpart_is_active(pi, e);
        const float hi = pi->h;
        const float hig2 = hi * hi * kernel_gamma2;
        const float pix[3] = {(float)(pi->x[0] - c->loc[0]),
            (float)(pi->x[1] - c->loc[1]),
            (float)(pi->x[2] - c->loc[2])};
        
        /* Get i particle time-step */
        const integertime_t ti_step = get_integer_timestep(pi->time_bin);
        const integertime_t ti_begin = get_integer_time_begin(e->ti_current - 1, pi->time_bin);
        double dti;
        if (with_cosmology) {
            dti = cosmology_get_delta_time(e->cosmology, ti_begin,
                                           ti_begin + ti_step);
        } else {
            dti = get_timestep(pi->time_bin, e->time_base);
        }
        
        
        /* Loop over the parts in cj. */
        for (int pjd = pid + 1; pjd < count; pjd++) {
            
            /* Get a pointer to the jth particle. */
            struct dmpart *restrict pj = &dmparts[pjd];
            
            /* Skip inhibited particles. */
            if (dmpart_is_inhibited(pj, e)) continue;
            
            /* Get j particle time-step */
            const integertime_t ti_step_j = get_integer_timestep(pj->time_bin);
            const integertime_t ti_begin_j = get_integer_time_begin(e->ti_current - 1, pj->time_bin);
            double dtj;
            if (with_cosmology) {
                dtj = cosmology_get_delta_time(e->cosmology, ti_begin_j,
                                               ti_begin_j + ti_step_j);
            } else {
                dtj = get_timestep(pj->time_bin, e->time_base);
            }
            
            const float hj = pj->h;
            const float hjg2 = hj * hj * kernel_gamma2;
            const int pj_active = dmpart_is_active(pj, e);
            
            /* Compute the pairwise distance. */
            const float pjx[3] = {(float)(pj->x[0] - c->loc[0]),
                (float)(pj->x[1] - c->loc[1]),
                (float)(pj->x[2] - c->loc[2])};
            float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
            const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            
            const int doi = pi_active && (r2 < hig2);
            const int doj = pj_active && (r2 < hjg2);
            
            /* Hit or miss? */
            if (doi && doj) {
                
                runner_iact_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, dti, dtj, ti_begin, sidm_props, us);
                
            } else if (doi) {
                
                runner_iact_nonsym_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, dti, dtj, ti_begin, sidm_props, us);
                
            } else if (doj) {
                
                dx[0] = -dx[0];
                dx[1] = -dx[1];
                dx[2] = -dx[2];
                
                runner_iact_nonsym_dark_matter_sidm(r2, dx, hj, hi, pj, pi, a, H, dtj, dti, ti_begin, sidm_props, us);
                
            }
        } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
    
    TIMER_TOC(TIMER_DOSELF);
}


/**
 * @brief
 *
 * @param r The thread #runner.
 * @param ci First #cell.
 * @param cj Second #cell.
 * @param timer Are we timing this?
 */
void runner_do_pair_dark_matter_sidm(struct runner *r, struct cell *ci,
                                struct cell *cj, int timer) {
    
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    const int with_cosmology = e->policy & engine_policy_cosmology;
    const struct unit_system *us = e->internal_units;
    const struct sidm_props *sidm_props = e->sidm_properties;
    
    TIMER_TIC;
    
    /* Anything to do here? */
    if (!cell_is_active_dark_matter(ci, e) && !cell_is_active_dark_matter(cj, e)) return;
    
    const int count_i = ci->dark_matter.count;
    const int count_j = cj->dark_matter.count;
    struct dmpart *restrict dmparts_i = ci->dark_matter.parts;
    struct dmpart *restrict dmparts_j = cj->dark_matter.parts;
    
    /* Cosmological terms */
    const float a = cosmo->a;
    const float H = cosmo->H;
    
    /* Get the relative distance between the pairs, wrapping. */
    double shift[3] = {0.0, 0.0, 0.0};
    for (int k = 0; k < 3; k++) {
        if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
            shift[k] = e->s->dim[k];
        else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
            shift[k] = -e->s->dim[k];
    }
    
    /* Loop over the parts in ci. */
    for (int pid = 0; pid < count_i; pid++) {
        
        /* Get a hold of the ith part in ci. */
        struct dmpart *restrict pi = &dmparts_i[pid];
        
        /* Skip inhibited particles. */
        if (dmpart_is_inhibited(pi, e)) continue;
        
        /* Get i particle time-step */
        const integertime_t ti_step = get_integer_timestep(pi->time_bin);
        const integertime_t ti_begin = get_integer_time_begin(e->ti_current - 1, pi->time_bin);
        double dti;
        if (with_cosmology) {
            dti = cosmology_get_delta_time(e->cosmology, ti_begin,
                                           ti_begin + ti_step);
        } else {
            dti = get_timestep(pi->time_bin, e->time_base);
        }
        
        const int pi_active = dmpart_is_active(pi, e);
        const float hi = pi->h;
        const float hig2 = hi * hi * kernel_gamma2;
        const float pix[3] = {(float)(pi->x[0] - (cj->loc[0] + shift[0])),
            (float)(pi->x[1] - (cj->loc[1] + shift[1])),
            (float)(pi->x[2] - (cj->loc[2] + shift[2]))};
        
        /* Loop over the parts in cj. */
        for (int pjd = 0; pjd < count_j; pjd++) {
            
            /* Get a pointer to the jth particle. */
            struct dmpart *restrict pj = &dmparts_j[pjd];
            
            /* Skip inhibited particles. */
            if (dmpart_is_inhibited(pj, e)) continue;
            
            /* Get j particle time-step */
            const integertime_t ti_step_j = get_integer_timestep(pj->time_bin);
            const integertime_t ti_begin_j = get_integer_time_begin(e->ti_current - 1, pj->time_bin);
            double dtj;
            if (with_cosmology) {
                dtj = cosmology_get_delta_time(e->cosmology, ti_begin_j,
                                               ti_begin_j + ti_step_j);
            } else {
                dtj = get_timestep(pj->time_bin, e->time_base);
            }
            
            const float hj = pj->h;
            const float hjg2 = hj * hj * kernel_gamma2;
            const int pj_active = dmpart_is_active(pj, e);
            
            /* Compute the pairwise distance. */
            const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                (float)(pj->x[1] - cj->loc[1]),
                (float)(pj->x[2] - cj->loc[2])};
            float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
            const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            
            /* Hit or miss? */
            if (r2 < hig2 && pi_active) {
                
                runner_iact_nonsym_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, dti, dtj, ti_begin, sidm_props, us);
                
            }
            if (r2 < hjg2 && pj_active) {
                
                dx[0] = -dx[0];
                dx[1] = -dx[1];
                dx[2] = -dx[2];
                
                runner_iact_nonsym_dark_matter_sidm(r2, dx, hj, hi, pj, pi, a, H, dtj, dti, ti_begin, sidm_props, us);
                
            }
        } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
    
    TIMER_TOC(TIMER_DOPAIR);


}


