/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2013 Pedro Gonnet (pedro.gonnet@durham.ac.uk)
 *               2016 Matthieu Schaller (matthieu.schaller@durham.ac.uk)
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
#include "../config.h"

/* This object's header. */
#include "runner_doiact_sidm.h"

/* Local includes. */
#include "active.h"
#include "cell.h"
#include "sidm.h"
#include "sidm_iact.h"
#include "inline.h"
#include "part.h"
#include "space_getsid.h"
#include "timers.h"


/**
 * @brief Computes the interaction of all the particles in a cell with all the
 * other ones.
 *
 * This function calculates the probability of DM-DM interactions
 *
 * @param r The #runner.
 * @param c The #cell we are working on.
 * @param timer Are we timing this ?
 */
void runner_doself_sidm(struct runner *r, struct cell *c, int timer) {

    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    
    TIMER_TIC;
    
    /* Cosmological terms */
    const float a = cosmo->a;
    const float H = cosmo->H;
        
    /* Anything to do here? */
    if (!cell_is_active_gravity(c, e)) return;

    /* Num parts */
    const int gcount = c->grav.count;
    struct part *restrict parts = c->grav.parts;

    /* Loop over all particles in ci... */
    for (int pid = 0; pid < gcount; pid++) {
        
        /* Get a pointer to the ith particle. */
        struct part *restrict pi = &parts[pid];
        
        /* Skip inactive particles */
        if (!part_is_active(pi, e)) continue;

        /* Loop over every other particle in the cell. */
        for (int pjd = 0; pjd < gcount_padded; pjd++) {
            
            /* No self interaction */
            if (pid == pjd) continue;
            
            /* Get a pointer to the jth particle. */
            struct part *restrict pj = &parts[pjd];

            /* Hit or miss? */
            if (r2 < hig2) {
                runner_iact_sidm(r2, dx, hi, hj, pi, pj, a, H);
            }

        } /* loop over the parts in ci. */
    } /* loop over the other parts in ci. */
    
    if (timer) TIMER_TOC(timer_doself_sidm);
}

/**
 * @brief Computes the interaction of all the particles in a cell with all the
 * particles of another cell.
 *
 * This function calculates the probability of DM-DM interactions.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The other #cell.
 */
void runner_dopair_sidm(struct runner *r, struct cell *ci, struct cell *cj) {
    
    const struct engine *restrict e = r->e;
    
    /* Get some other useful values. */
    const double hi_max = ci->hydro.h_max;
    const double hj_max = cj->hydro.h_max;
    const int count_i = ci->hydro.count;
    const int count_j = cj->hydro.count;
    struct part *restrict parts_i = ci->hydro.parts;
    struct part *restrict parts_j = cj->hydro.parts;

    
    /* Anything to do here? */
    if (count_i == 0 || count_j == 0) return;
    
    /* Anything to do here? */
    if (!cell_is_active_gravity(ci, e) && !cell_is_active_gravity(cj, e)) return;
    
    for (int pid = 0; pid < count_i; pid++) {
        
        struct part *pi = &parts_i[pid];
        
        /* Skip inactive particles */
        if (!part_is_active(pi, e)) continue;
        
        for (int pjd = 0; pjd < count_j; pjd++) {
            
            struct part *pj = &parts_i[pjd];
            
            /* Skip inactive particles */
            if (!part_is_active(pj, e)) continue;
            
            /* Hit or miss? */
            if (r2 < hig2) {
                runner_iact_sidm(r2, dx, hi, hj, pi, pj, a, H);
            }
        }
    }
    
    TIMER_TOC(timer_dopair_sidm);
}
