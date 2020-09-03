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
   runner_dopair_FUNCTION, runner_doself_FUNCTION,
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
void DOPAIR_DM(struct runner *r, struct cell *restrict ci,
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
            if (part_is_inhibited(pj, e)) continue;

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
 * @brief Compute the interactions within a cell (non-symmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param c The #cell.
 */
void DOSELF_DM(struct runner *r, struct cell *restrict c) {
    
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
        
        const int pi_active = dmpart_is_active(pi, e);
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
 * @brief Compute grouped sub-cell interactions for pairs
 *
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param parts_i The #part to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param count The number of particles in @c ind.
 * @param cj The second #cell.
 */
void DOSUB_PAIR_DM(struct runner *r, struct cell *restrict ci, struct cell *restrict cj, int count) {
    
    struct space *s = r->e->s;
    const struct engine *e = r->e;
    
    TIMER_TIC;
    
    /* Should we even bother? */
    if (!cell_is_active_dark_matter(ci, e) && !cell_is_active_dark_matter(cj, e)) return;
    if (ci->dark_matter.count == 0 || cj->dark_matter.count == 0) return;
    
    /* Get the type of pair and flip ci/cj if needed. */
    double shift[3];
    const int sid = space_getsid(s, &ci, &cj, shift);
    
    /* Recurse? */
    if (cell_can_recurse_in_pair_dark_matter_task(ci) &&
        cell_can_recurse_in_pair_dark_matter_task(cj)) {
        struct cell_split_pair *csp = &cell_split_pairs[sid];
        for (int k = 0; k < csp->count; k++) {
            const int pid = csp->pairs[k].pid;
            const int pjd = csp->pairs[k].pjd;
            if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL)
                DOSUB_PAIR_DM(r, ci->progeny[pid], cj->progeny[pjd], 0);
        }
    }
    
    /* Otherwise, compute the pair directly. */
    else if (cell_is_active_dark_matter(ci, e) || cell_is_active_dark_matter(cj, e)) {
        
        /* Make sure both cells are drifted to the current timestep. */
        if (!cell_are_dmpart_drifted(ci, e) || !cell_are_dmpart_drifted(cj, e))
            error("Interacting undrifted cells.");
        
        /* Compute the interactions. */
        DOPAIR_DM(r, ci, cj);
    }
    
    if (gettimer) TIMER_TOC(TIMER_DOSUB_PAIR_DM);
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
void DOSUB_SELF_DM(struct runner *r, struct cell *restrict ci, int count) {
    
    TIMER_TIC;
    
    /* Should we even bother? */
    if (ci->dark_matter.count == 0 || !cell_is_active_dark_matter(ci, r->e)) return;
    
    /* Recurse? */
    if (cell_can_recurse_in_self_dark_matter_task(ci)) {
        
        /* Loop over all progeny. */
        for (int k = 0; k < 8; k++)
            if (ci->progeny[k] != NULL) {
                DOSUB_SELF_DM(r, ci->progeny[k], 0);
                for (int j = k + 1; j < 8; j++)
                    if (ci->progeny[j] != NULL)
                        DOSUB_PAIR_DM(r, ci->progeny[k], ci->progeny[j], 0);
            }
    }
    
    /* Otherwise, compute self-interaction. */
    else {
        
        /* Drift the cell to the current timestep if needed. */
        if (!cell_are_dmpart_drifted(ci, r->e)) error("Interacting undrifted cell.");
        
        DOSELF_DM(r, ci);
    }
    
    if (gettimer) TIMER_TOC(TIMER_DOSUB_SELF_DM);
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
void DOSELF_SUBSET_DM(struct runner *r, struct cell *restrict ci,
                   struct dmpart *restrict dmparts, int *restrict ind, int count) {
    
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    
    TIMER_TIC;
    
    /* Cosmological terms */
    const float a = cosmo->a;
    const float H = cosmo->H;
    
    const int count_i = ci->dark_matter.count;
    struct dmpart *restrict dmparts_j = ci->dark_matter.parts;
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
    
    TIMER_TOC(timer_doself_subset_dm);
}

/**
 * @brief Compute the interactions between a cell pair, but only for the
 *      given indices in ci.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param parts_i The #part to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param count The number of particles in @c ind.
 * @param cj The second #cell.
 * @param sid The direction of the pair.
 * @param flipped Flag to check whether the cells have been flipped or not.
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR_SUBSET_DM(struct runner *r, struct cell *restrict ci,
                   struct part *restrict parts_i, int *restrict ind, int count,
                   struct cell *restrict cj, const int sid, const int flipped,
                   const double *shift) {
    
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    
    TIMER_TIC;
    
    const int count_j = cj->dark_matter.count;
    struct dmpart *restrict dmparts_j = cj->dark_matter.parts;
    
    /* Cosmological terms */
    const float a = cosmo->a;
    const float H = cosmo->H;
    
    /* Pick-out the sorted lists. */
    const struct sort_entry *sort_j = cell_get_dark_matter_sorts(cj, sid);
    const float dxj = cj->dark_matter.dx_max_sort;
    
    /* Parts are on the left? */
    if (!flipped) {
        
        /* Loop over the parts_i. */
        for (int pid = 0; pid < count; pid++) {
            
            /* Get a hold of the ith part in ci. */
            struct dmpart *restrict pi = &dmparts_i[ind[pid]];
            const double pix = pi->x[0] - (shift[0]);
            const double piy = pi->x[1] - (shift[1]);
            const double piz = pi->x[2] - (shift[2]);
            const float hi = pi->h;
            const float hig2 = hi * hi * kernel_gamma2;
            const double di = hi * kernel_gamma + dxj + pix * runner_shift[sid][0] +
            piy * runner_shift[sid][1] + piz * runner_shift[sid][2];
            
            /* Loop over the parts in cj. */
            for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {
                
                /* Get a pointer to the jth particle. */
                struct dmpart *restrict pj = &dmparts_j[sort_j[pjd].i];
                
                /* Skip inhibited particles. */
                if (dmpart_is_inhibited(pj, e)) continue;
                
                const float hj = pj->h;
                const double pjx = pj->x[0];
                const double pjy = pj->x[1];
                const double pjz = pj->x[2];
                
                /* Compute the pairwise distance. */
                float dx[3] = {(float)(pix - pjx), (float)(piy - pjy),
                    (float)(piz - pjz)};
                const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
                
                /* Hit or miss? */
                if (r2 < hig2) {
                    
                    IACT_NONSYM_DM(r2, dx, hi, hj, pi, pj, a, H);
                }
            } /* loop over the parts in cj. */
        }   /* loop over the parts in ci. */
    }
    
    /* Parts are on the right. */
    else {
        
        /* Loop over the parts_i. */
        for (int pid = 0; pid < count; pid++) {
            
            /* Get a hold of the ith part in ci. */
            struct dmpart *restrict pi = &dmparts_i[ind[pid]];
            const double pix = pi->x[0] - (shift[0]);
            const double piy = pi->x[1] - (shift[1]);
            const double piz = pi->x[2] - (shift[2]);
            const float hi = pi->h;
            const float hig2 = hi * hi * kernel_gamma2;
            const double di = -hi * kernel_gamma - dxj + pix * runner_shift[sid][0] +
            piy * runner_shift[sid][1] + piz * runner_shift[sid][2];
            
            /* Loop over the parts in cj. */
            for (int pjd = count_j - 1; pjd >= 0 && di < sort_j[pjd].d; pjd--) {
                
                /* Get a pointer to the jth particle. */
                struct dmpart *restrict pj = &dmparts_j[sort_j[pjd].i];
                
                /* Skip inhibited particles. */
                if (dmpart_is_inhibited(pj, e)) continue;
                
                const float hj = pj->h;
                const double pjx = pj->x[0];
                const double pjy = pj->x[1];
                const double pjz = pj->x[2];
                
                /* Compute the pairwise distance. */
                float dx[3] = {(float)(pix - pjx), (float)(piy - pjy),
                    (float)(piz - pjz)};
                const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
         
                /* Hit or miss? */
                if (r2 < hig2) {
                    
                    IACT_NONSYM_DM(r2, dx, hi, hj, pi, pj, a, H);

                }
            } /* loop over the parts in cj. */
        }   /* loop over the parts in ci. */
    }
    
    TIMER_TOC(timer_dopair_subset);
}

