/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2020 Camila Correa (camila.correa@uva.nl)
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
#include "inline.h"
#include "part.h"
#include "space_getsid.h"
#include "timers.h"

/* sidm inclues. */
#include "sidm.h"
#include "sidm_iact.h"
#include "sidm_properties.h"


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
void runner_doself_sidm(struct runner *r, struct cell *c) {

    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    const int with_cosmology = (e->policy & engine_policy_cosmology);
    const struct gravity_props *grav_props = e->gravity_properties;
    const struct unit_system *us = e->internal_units;
    const struct sidm_props *sidm_props = e->sidm_properties;

    
    /* Cosmological terms */
    const float a = cosmo->a;
    const float H = cosmo->H;
        
    /* Anything to do here? */
    if (!cell_is_active_gravity(c, e)) return;

    /* Num parts in cell */
    const int gcount = c->grav.count;
    struct gpart *gparts = c->grav.parts;

    /* Loop over all particles in cell... */
    for (int pid = 0; pid < gcount; pid++) {
        
        /* Get a pointer to the ith particle. */
        struct gpart *gpi = &gparts[pid];
        
        /* Skip inactive particles */
        if (!gpart_is_active(gpi, e)) continue;
        
        const float hi = gravity_get_softening(gpi, grav_props);
        
        /* Loop over every other particle in the cell. */
        for (int pjd = 0; pjd < gcount; pjd++) {
            
            /* No self interaction */
            if (pid == pjd) continue;
            
            /* Get a pointer to the jth particle. */
            struct gpart *gpj = &gparts[pjd];
            
            const integertime_t ti_step = get_integer_timestep(gpj->time_bin);
            const integertime_t ti_begin = get_integer_time_begin(e->ti_current - 1, gpj->time_bin);
            
            /* Get particle time-step */
            double dt;
            if (with_cosmology) {
                dt = cosmology_get_delta_time(e->cosmology, ti_begin,
                                                   ti_begin + ti_step);
            } else {
                dt = get_timestep(gpj->time_bin, e->time_base);
            }
            
            const double Gyr_in_cgs = 1e9 * 365.25 * 24. * 3600.;
            const double time_to_cgs = units_cgs_conversion_factor(us, UNIT_CONV_TIME);
            const double conversion_factor = time_to_cgs / Gyr_in_cgs;
            const double dt_Gyr = dt * conversion_factor;
            
            const float hj = gravity_get_softening(gpj, grav_props);
            const float hj2 = hj * hj;
            
            /* Compute the pairwise distance. */
            float dx[3] = {gpj->x[0] - gpi->x[0], gpj->x[1] - gpi->x[1], gpj->x[2] - gpi->x[2]};
            const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

            /* Hit or miss? */
            if (r2 < hj2) {
                runner_iact_sidm(r2, dx, hi, hj, gpi, gpj, a, H, dt_Gyr, sidm_props);
            }

        } /* loop over the parts in cell. */
    } /* loop over the other parts in cell. */
    
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
    
}



