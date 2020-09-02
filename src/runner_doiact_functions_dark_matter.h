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

/* Before including this file, define FUNCTION, which is the
   name of the interaction function. This creates the interaction functions
   runner_dopair_FUNCTION, runner_dopair_FUNCTION_naive, runner_doself_FUNCTION,
   and runner_dosub_FUNCTION calling the pairwise interaction function
   runner_iact_FUNCTION. */

#include "runner_doiact_dark_matter.h"

/**
 * @brief Compute the interactions between a cell pair (non-symmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
void DOPAIR1_DM_NAIVE(struct runner *r, struct cell *restrict ci,
                   struct cell *restrict cj) {
    
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    
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
        
        const int pi_active = dark_matter_is_active(pi, e);
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
            if (part_is_inhibited(pj, e)) continue;
            
            const float hj = pj->h;
            const float hjg2 = hj * hj * kernel_gamma2;
            const int pj_active = dark_matter_is_active(pj, e);
            
            /* Compute the pairwise distance. */
            const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                (float)(pj->x[1] - cj->loc[1]),
                (float)(pj->x[2] - cj->loc[2])};
            float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
            const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            
            /* Hit or miss? */
            if (r2 < hig2 && pi_active) {
                
                IACT_NONSYM_DM(r2, dx, hi, hj, pi, pj, a, H);
            }
            if (r2 < hjg2 && pj_active) {
                
                dx[0] = -dx[0];
                dx[1] = -dx[1];
                dx[2] = -dx[2];
                
                IACT_NONSYM_DM(r2, dx, hj, hi, pj, pi, a, H);
            }
        } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
    
    TIMER_TOC(TIMER_DOPAIR);
}

/**
 * @brief Compute the interactions between a cell pair (symmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
void DOPAIR2_DM_NAIVE(struct runner *r, struct cell *restrict ci,
                   struct cell *restrict cj) {
    
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    
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
            
            const int pj_active = dmpart_is_active(pj, e);
            const float hj = pj->h;
            const float hjg2 = hj * hj * kernel_gamma2;
            
            /* Compute the pairwise distance. */
            const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                (float)(pj->x[1] - cj->loc[1]),
                (float)(pj->x[2] - cj->loc[2])};
            float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
            const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            
            /* Hit or miss? */
            if (r2 < hig2 || r2 < hjg2) {
                
                if (pi_active && pj_active) {
                    
                    IACT_DM(r2, dx, hi, hj, pi, pj, a, H);
                    
                } else if (pi_active) {
                    
                    IACT_NONSYM_DM(r2, dx, hi, hj, pi, pj, a, H);
                    
                } else if (pj_active) {
                    
                    dx[0] = -dx[0];
                    dx[1] = -dx[1];
                    dx[2] = -dx[2];
                    
                    IACT_NONSYM_DM(r2, dx, hj, hi, pj, pi, a, H);
                }
            }
        } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
    
    TIMER_TOC(TIMER_DOPAIR);
}

/**
 * @brief Compute the interactions within a cell (non-symmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param c The #cell.
 */
void DOSELF1_DM_NAIVE(struct runner *r, struct cell *restrict c) {
    
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    
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
        
        const int pi_active = dark_matter_is_active(pi, e);
        const float hi = pi->h;
        const float hig2 = hi * hi * kernel_gamma2;
        const float pix[3] = {(float)(pi->x[0] - c->loc[0]),
            (float)(pi->x[1] - c->loc[1]),
            (float)(pi->x[2] - c->loc[2])};
        
        /* Loop over the parts in cj. */
        for (int pjd = pid + 1; pjd < count; pjd++) {
            
            /* Get a pointer to the jth particle. */
            struct part *restrict pj = &dmparts[pjd];
            
            /* Skip inhibited particles. */
            if (dmpart_is_inhibited(pj, e)) continue;
            
            const float hj = pj->h;
            const float hjg2 = hj * hj * kernel_gamma2;
            const int pj_active = dark_matter_is_active(pj, e);
            
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
                
                IACT_DM(r2, dx, hi, hj, pi, pj, a, H);

            } else if (doi) {
                
                IACT_NONSYM_DM(r2, dx, hi, hj, pi, pj, a, H);

            } else if (doj) {
                
                dx[0] = -dx[0];
                dx[1] = -dx[1];
                dx[2] = -dx[2];
                
                IACT_NONSYM_DM(r2, dx, hj, hi, pj, pi, a, H);

            }
        } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
    
    TIMER_TOC(TIMER_DOSELF);
}

/**
 * @brief Compute the interactions within a cell (symmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param c The #cell.
 */
void DOSELF2_DM_NAIVE(struct runner *r, struct cell *restrict c) {
    
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    
    TIMER_TIC;
    
    /* Anything to do here? */
    if (!cell_is_active_dark_matter(c, e)) return;
    
    /* Cosmological terms */
    const float a = cosmo->a;
    const float H = cosmo->H;
    
    const int count = c->dark_matter.count;
    struct part *restrict parts = c->dark_matter.parts;
    
    /* Loop over the parts in ci. */
    for (int pid = 0; pid < count; pid++) {
        
        /* Get a hold of the ith part in ci. */
        struct dmpart *restrict pi = &dmparts[pid];
        
        /* Skip inhibited particles. */
        if (dmpart_is_inhibited(pi, e)) continue;
        
        const int pi_active = dark_matter_is_active(pi, e);
        const float hi = pi->h;
        const float hig2 = hi * hi * kernel_gamma2;
        const float pix[3] = {(float)(pi->x[0] - c->loc[0]),
            (float)(pi->x[1] - c->loc[1]),
            (float)(pi->x[2] - c->loc[2])};
        
        /* Loop over the parts in cj. */
        for (int pjd = pid + 1; pjd < count; pjd++) {
            
            /* Get a pointer to the jth particle. */
            struct dmpart *restrict pj = &dmparts[pjd];
            
            /* Skip inhibited particles. */
            if (dmpart_is_inhibited(pj, e)) continue;
            
            const float hj = pj->h;
            const float hjg2 = hj * hj * kernel_gamma2;
            const int pj_active = dark_matter_is_active(pj, e);
            
            /* Compute the pairwise distance. */
            const float pjx[3] = {(float)(pj->x[0] - c->loc[0]),
                (float)(pj->x[1] - c->loc[1]),
                (float)(pj->x[2] - c->loc[2])};
            float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
            const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            
            const int doi = pi_active && ((r2 < hig2) || (r2 < hjg2));
            const int doj = pj_active && ((r2 < hig2) || (r2 < hjg2));
            
            /* Hit or miss? */
            if (doi && doj) {
                
                IACT_DM(r2, dx, hi, hj, pi, pj, a, H);

            } else if (doi) {
                
                IACT_NONSYM_DM(r2, dx, hi, hj, pi, pj, a, H);

            } else if (doj) {
                
                dx[0] = -dx[0];
                dx[1] = -dx[1];
                dx[2] = -dx[2];
                
                IACT_NONSYM_DM(r2, dx, hj, hi, pj, pi, a, H);

            }
        } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
    
    TIMER_TOC(TIMER_DOSELF);
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * Version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param parts_i The #part to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param count The number of particles in @c ind.
 * @param cj The second #cell.
 */
void DOPAIR_SUBSET_DM_NAIVE(struct runner *r, struct cell *restrict ci,
                         struct dmpart *restrict dmparts_i, int *restrict ind,
                         int count, struct cell *restrict cj) {
    
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    
    TIMER_TIC;
    
    const int count_j = cj->dark_matter.count;
    struct dmpart *restrict dmparts_j = cj->dark_matter.parts;
    
    /* Cosmological terms */
    const float a = cosmo->a;
    const float H = cosmo->H;
    
    /* Loop over the parts_i. */
    for (int pid = 0; pid < count; pid++) {
        
        /* Get a hold of the ith part in ci. */
        struct dmpart *restrict pi = &dmparts_i[ind[pid]];
        double pix[3];
        for (int k = 0; k < 3; k++) pix[k] = pi->x[k] - shift[k];
        const float hi = pi->h;
        const float hig2 = hi * hi * kernel_gamma2;
        
        /* Loop over the parts in cj. */
        for (int pjd = 0; pjd < count_j; pjd++) {
            
            /* Get a pointer to the jth particle. */
            struct dmpart *restrict pj = &dmparts_j[pjd];
            
            /* Skip inhibited particles. */
            if (dmpart_is_inhibited(pj, e)) continue;
            
            /* Compute the pairwise distance. */
            float r2 = 0.0f;
            float dx[3];
            for (int k = 0; k < 3; k++) {
                dx[k] = pix[k] - pj->x[k];
                r2 += dx[k] * dx[k];
            }
            
            /* Hit or miss? */
            if (r2 < hig2) {
                
                IACT_NONSYM_DM(r2, dx, hi, pj->h, pi, pj, a, H);
                
            }
        } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
    
    TIMER_TOC(timer_dopair_subset_dm_naive);
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param parts The #part to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param count The number of particles in @c ind.
 */
void DOSELF_SUBSET_DM_NAIVE(struct runner *r, struct cell *restrict ci,
                   struct dmpart *restrict dmparts, int *restrict ind, int count) {
    
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    
    TIMER_TIC;
    
    /* Cosmological terms */
    const float a = cosmo->a;
    const float H = cosmo->H;
    
    const int count_i = ci->dark_matter.count;
    struct dmpart *restrict parts_j = ci->dark_matter.parts;
    /* Loop over the parts in ci. */
    for (int pid = 0; pid < count; pid++) {
        
        /* Get a hold of the ith part in ci. */
        struct dmpart *pi = &dmparts[ind[pid]];
        const float pix[3] = {(float)(pi->x[0] - ci->loc[0]),
            (float)(pi->x[1] - ci->loc[1]),
            (float)(pi->x[2] - ci->loc[2])};
        const float hi = pi->h;
        const float hig2 = hi * hi * kernel_gamma2;
        
        /* Loop over the parts in cj. */
        for (int pjd = 0; pjd < count_i; pjd++) {
            
            /* Get a pointer to the jth particle. */
            struct dmpart *restrict pj = &dmparts_j[pjd];
            
            /* Skip inhibited particles. */
            if (dmpart_is_inhibited(pj, e)) continue;
            
            const float hj = pj->h;
            
            /* Compute the pairwise distance. */
            const float pjx[3] = {(float)(pj->x[0] - ci->loc[0]),
                (float)(pj->x[1] - ci->loc[1]),
                (float)(pj->x[2] - ci->loc[2])};
            float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
            const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            
            /* Hit or miss? */
            if (r2 > 0.f && r2 < hig2) {
                
                IACT_NONSYM_DM(r2, dx, hi, hj, pi, pj, a, H);
            }
        } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
    
    TIMER_TOC(timer_doself_subset);
}







