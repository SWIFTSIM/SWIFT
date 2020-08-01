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
 * @brief Density calculation for dark matter particles between two cells.
 *
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
void runner_dopair_dark_matter_density(struct runner *r, struct cell *restrict ci,
                   struct cell *restrict cj) {
    
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    
    /* Anything to do here? */
    if (!cell_is_active_gravity(ci, e)) return;
    if (ci->grav.count == 0) return;
    
    /* Num parts in cells */
    const int gcount_i = ci->grav.count;
    const int gcount_j = cj->grav.count;
    
    /* Select part in ci. */
    struct gpart *gparts_i = ci->grav.parts;
    struct gpart *gparts_j = cj->grav.parts;

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
        struct gpart *gpi = &gparts_i[pid];
        
        /* We care only about dark matter particles */
        if (gpi->type != swift_type_dark_matter) continue;
        
        const int gpi_active = gpart_is_active(gpi, e);
        const float hi = gpi->sidm_data.h_sidm;
        const float hig2 = hi * hi * kernel_gamma2;
        const float gpix[3] = {(float)(gpi->x[0] - (cj->loc[0] + shift[0])),
            (float)(gpi->x[1] - (cj->loc[1] + shift[1])),
            (float)(gpi->x[2] - (cj->loc[2] + shift[2]))};
        
        /* Loop over the parts in cj. */
        for (int pjd = 0; pjd < count_j; pjd++) {
            
            /* Get a pointer to the jth particle. */
            struct gpart *restrict gpj = &gparts_j[pjd];
            
            /* Use only dark matter particles */
            if (gpj->type != swift_type_dark_matter) continue;
            
            const float hj = gpj->sidm_data.h_sidm;
            const float hjg2 = hj * hj * kernel_gamma2;
            const int gpj_active = gpart_is_active(gpj, e);
            
            /* Compute the pairwise distance. */
            const float gpjx[3] = {(float)(gpj->x[0] - cj->loc[0]),
                (float)(gpj->x[1] - cj->loc[1]),
                (float)(gpj->x[2] - cj->loc[2])};
            float dx[3] = {gpix[0] - gpjx[0], gpix[1] - gpjx[1], gpix[2] - gpjx[2]};
            const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            
            /* Hit or miss? */
            if (r2 < hig2 && pi_active) {
                
                runner_iact_nonsym_dm_density(r2, dx, hi, hj, pi, pj, a, H);

            }
            if (r2 < hjg2 && pj_active) {
                
                dx[0] = -dx[0];
                dx[1] = -dx[1];
                dx[2] = -dx[2];
                
                runner_iact_nonsym_dm_density(r2, dx, hj, hi, pj, pi, a, H);
            }
        } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
}

/**
 * @brief Compute the interactions within a cell (non-symmetric case).
 *
 * @param r The #runner.
 * @param c The #cell.
 */
void runner_doself_dark_matter_density(struct runner *r, struct cell *restrict c) {
    
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    
    /* Anything to do here? */
    if (!cell_is_active_gravity(ci, e)) return;
    if (ci->grav.count == 0) return;
    
    /* Cosmological terms */
    const float a = cosmo->a;
    const float H = cosmo->H;
    
    /* Num parts in cell */
    const int gcount = c->grav.count;
    struct gpart *gparts = c->grav.parts;
    
    /* Loop over the parts in ci. */
    for (int pid = 0; pid < gcount; pid++) {
        
        /* Get a hold of the ith part in ci. */
        struct gpart *restrict gpi = &gparts[pid];
        
        const int gpi_active = gpart_is_active(gpi, e);
        const float hi = gpi->sidm_data.h_sidm;
        const float hig2 = hi * hi * kernel_gamma2;
        const float gpix[3] = {(float)(gpi->x[0] - c->loc[0]),
            (float)(gpi->x[1] - c->loc[1]),
            (float)(gpi->x[2] - c->loc[2])};
        
        /* Loop over the parts in cj. */
        for (int pjd = pid + 1; pjd < count; pjd++) {
            
            /* Get a pointer to the jth particle. */
            struct gpart *restrict gpj = &gparts[pjd];
            
            const float hj = gpj->sidm_data.h_sidm;
            const float hjg2 = hj * hj * kernel_gamma2;
            const int gpj_active = gpart_is_active(gpj, e);
            
            /* Compute the pairwise distance. */
            const float gpjx[3] = {(float)(gpj->x[0] - c->loc[0]),
                (float)(gpj->x[1] - c->loc[1]),
                (float)(gpj->x[2] - c->loc[2])};
            float dx[3] = {gpix[0] - gpjx[0], gpix[1] - gpjx[1], gpix[2] - gpjx[2]};
            const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            
            const int doi = gpi_active && (r2 < hig2);
            const int doj = gpj_active && (r2 < hjg2);
            
            /* Hit or miss? */
            if (doi && doj) {
                
                runner_iact_dark_matter_density(r2, dx, hi, hj, pi, pj, a, H);

            } else if (doi) {
                
                runner_iact_nonsym_dark_matter_density(r2, dx, hi, hj, pi, pj, a, H);

            } else if (doj) {
                
                dx[0] = -dx[0];
                dx[1] = -dx[1];
                dx[2] = -dx[2];
                
                runner_iact_nonsym_dark_matter_density(r2, dx, hj, hi, pj, pi, a, H);

            }
        } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
}

/**
 * @brief Computes the interaction of all the particles within a cell
 *
 * This function calculates the probability of DM-DM interactions
 *
 * @param r The #runner.
 * @param c The #cell we are working on.
 */
void runner_doself_sidm(struct runner *r, struct cell *c) {

    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    const int with_cosmology = (e->policy & engine_policy_cosmology);
    const struct unit_system *us = e->internal_units;
    const struct sidm_props *sidm_props = e->sidm_properties;

    /* Cosmological terms */
    const float a = cosmo->a;
    const float H = cosmo->H;
        
    /* Anything to do here? */
    if (!cell_is_active_gravity(c, e)) return;
    if (c->grav.count == 0) return;
    
    /* Num parts in cell */
    const int gcount = c->grav.count;
    struct gpart *gparts = c->grav.parts;

    /* Loop over all particles in cell... */
    for (int pid = 0; pid < gcount; pid++) {
        
        /* Get a pointer to the ith particle. */
        struct gpart *gpi = &gparts[pid];
        
        /* We care only about dark matter particles here */
        if (gpi->type != swift_type_dark_matter) continue;

        /* Skip inactive particles */
        if (!gpart_is_active(gpi, e)) continue;
        
        /* Get particle time-step */
        const integertime_t ti_step = get_integer_timestep(gpi->time_bin);
        const integertime_t ti_begin = get_integer_time_begin(e->ti_current - 1, gpi->time_bin);
        double dt;
        if (with_cosmology) {
            dt = cosmology_get_delta_time(e->cosmology, ti_begin,
                                          ti_begin + ti_step);
        } else {
            dt = get_timestep(gpi->time_bin, e->time_base);
        }
        
        /* Loop over every other particle in the cell. */
        for (int pjd = pid + 1; pjd < gcount; pjd++) {
            
            /* No self interaction */
            if (pid == pjd) continue;
            
            /* Get a pointer to the jth particle. */
            struct gpart *gpj = &gparts[pjd];
            
            /* We care only about dark matter particles ... */
            if (gpj->type != swift_type_dark_matter) continue;
            
            /* Compute the pairwise distance. */
            float dx[3] = {gpi->x[0] - gpj->x[0], gpi->x[1] - gpj->x[1], gpi->x[2] - gpj->x[2]};
            const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            const float h_SI = gpj->sidm_data.h_sidm;
            const float h_SI2 = h_SI * h_SI;
            
            /* Hit or miss? */
            if (r2 < h_SI2) {
                runner_iact_sidm(h_SI, gpi, gpj, a, H, dt, ti_begin, sidm_props, us);
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
    
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    const int with_cosmology = (e->policy & engine_policy_cosmology);
    const struct unit_system *us = e->internal_units;
    const struct sidm_props *sidm_props = e->sidm_properties;
    
    /* Cosmological terms */
    const float a = cosmo->a;
    const float H = cosmo->H;
    
    /* Anything to do here? */
    if (!cell_is_active_gravity(ci, e)) return;
    if (ci->grav.count == 0) return;
    
    /* Num parts in cells */
    const int gcount_i = ci->grav.count;
    const int gcount_j = cj->grav.count;
    
    /* Select part in ci. */
    struct gpart *gparts_i = ci->grav.parts;
    struct gpart *gparts_j = cj->grav.parts;

    /* Loop over all particles in cell... */
    for (int pid = 0; pid < gcount_i; pid++) {
        
        /* Get a hold of the ith part in ci. */
        struct gpart *gpi = &gparts_i[pid];
        
        /* We care only about dark matter particles */
        if (gpi->type != swift_type_dark_matter) continue;
        
        /* Skip inactive particles */
        if (!gpart_is_active(gpi, e)) continue;
        
        /* Get particle time-step */
        const integertime_t ti_step = get_integer_timestep(gpi->time_bin);
        const integertime_t ti_begin = get_integer_time_begin(e->ti_current - 1, gpi->time_bin);
        double dt;
        if (with_cosmology) {
            dt = cosmology_get_delta_time(e->cosmology, ti_begin,
                                          ti_begin + ti_step);
        } else {
            dt = get_timestep(gpi->time_bin, e->time_base);
        }
        
        /* Loop over every other particle in the cell. */
        for (int pjd = 0; pjd < gcount_j; pjd++) {
            
            /* No self interaction */
            if (pid == pjd) continue;
            
            /* Get a pointer to the jth particle. */
            struct gpart *gpj = &gparts_j[pjd];
            
            /* We care only about dark matter particles */
            if (gpj->type != swift_type_dark_matter) continue;
            
            /* Compute the pairwise distance. */
            float dx[3] = {gpi->x[0] - gpj->x[0], gpi->x[1] - gpj->x[1], gpi->x[2] - gpj->x[2]};
            const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            const float h_SI = gpj->sidm_data.h_sidm;
            const float h_SI2 = h_SI * h_SI;
            
            /* Hit or miss? */
            if (r2 < h_SI2) {
                runner_iact_sidm(h_SI, gpi, gpj, a, H, dt, ti_begin, sidm_props, us);
            }
            
        } /* loop over the parts in cell cj. */
    } /* loop over the other parts in cell ci. */
    
}

/**
 * @brief Performs the kicks in momentum space from DM-DM interactions on all the active particles in a cell.
 *
 * @param r The runner thread.
 * @param c The cell.
 */
void runner_do_sidm_kick(struct runner *r, struct cell *c) {
    
    /* Anything to do here? */
    if (c->grav.count == 0) return;

    const struct engine *e = r->e;
    struct gpart *restrict gparts = c->grav.parts;
    const int gcount = c->grav.count;
    
    /* Loop over the gparts in this cell. */
    for (int k = 0; k < gcount; k++) {
        
        /* Get a handle on the part. */
        struct gpart *restrict gp = &gparts[k];
        
        /* Do it only to dark matter particles */
        if (gp->type != swift_type_dark_matter) continue;
        
        /* do the kick */
        communicate_sidm_kick_to_gpart(gp);
    }
}



