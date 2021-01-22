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

#include "runner_doiact_dark_matter.h"

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
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR_SUBSET_NAIVE(struct runner *r, struct cell *restrict ci,
                         struct dmpart *restrict dmparts_i, int *restrict ind,
                         int count, struct cell *restrict cj,
                         const double *shift) {
    
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
        const float hig2 = hi * hi * dm_kernel_gamma2;
        
        /* Loop over the parts in cj. */
        for (int pjd = 0; pjd < count_j; pjd++) {
            
            /* Get a pointer to the jth particle. */
            struct dmpart *restrict pj = &dmparts_j[pjd];
            
            /* Skip inhibited particles. */
            if (dmpart_is_inhibited(pj, e)) continue;
            
#ifdef SWIFT_DEBUG_CHECKS
            /* Check that particles have been drifted to the current time */
            if (pj->ti_drift != e->ti_current)
                error("DM Particle pj not drifted to current time");
#endif
            
            /* Compute the pairwise distance. */
            float r2 = 0.0f;
            float dx[3];
            for (int k = 0; k < 3; k++) {
                dx[k] = pix[k] - pj->x[k];
                r2 += dx[k] * dx[k];
            }
            
            /* Hit or miss? */
            if (r2 < hig2) {
                
                runner_iact_nonsym_dark_matter_density(r2, dx, hi, pj->h, pi, pj, a, H);
            }
        } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
    
    TIMER_TOC(timer_dopair_subset_naive);
}


/**
 * @brief Determine which version of DOPAIR_SUBSET needs to be called depending
 * on the
 * orientation of the cells or whether DOPAIR_SUBSET needs to be called at all.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param parts_i The #part to interact with @c cj.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param count The number of particles in @c ind.
 * @param cj The second #cell.
 */
void DOPAIR_SUBSET_BRANCH(struct runner *r, struct cell *restrict ci,
                          struct dmpart *restrict dmparts_i, int *restrict ind,
                          struct cell *restrict cj, int count) {
    
    const struct engine *e = r->e;
    
    /* Anything to do here? */
    if (cj->dark_matter.count == 0) return;
    
    /* Get the relative distance between the pairs, wrapping. */
    double shift[3] = {0.0, 0.0, 0.0};
    for (int k = 0; k < 3; k++) {
        if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
        shift[k] = e->s->dim[k];
        else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
        shift[k] = -e->s->dim[k];
    }
    
    DOPAIR_SUBSET_NAIVE(r, ci, dmparts_i, ind, count, cj, shift);
}


void DOSUB_SUBSET(struct runner *r, struct cell *ci, struct dmpart *dmparts,
                  int *ind, int count, struct cell *cj) {
    
    const struct engine *e = r->e;
    struct space *s = e->s;
    
    TIMER_TIC;
    
    /* Should we even bother? */
    if (!cell_is_active_dark_matter(ci, e) &&
        (cj == NULL || !cell_is_active_dark_matter(cj, e)))
    return;
    if (ci->dark_matter.count == 0 || (cj != NULL && cj->dark_matter.count == 0)) return;
    
    /* Find out in which sub-cell of ci the parts are. */
    struct cell *sub = NULL;
    if (ci->split) {
        for (int k = 0; k < 8; k++) {
            if (ci->progeny[k] != NULL) {
                if (&dmparts[ind[0]] >= &ci->progeny[k]->dark_matter.parts[0] &&
                    &dmparts[ind[0]] <
                    &ci->progeny[k]->dark_matter.parts[ci->progeny[k]->dark_matter.count]) {
                    sub = ci->progeny[k];
                    break;
                }
            }
        }
    }
    
    /* Is this a single cell? */
    if (cj == NULL) {
        
        /* Recurse? */
        if (cell_can_recurse_in_self_dark_matter_task(ci)) {
            
            /* Loop over all progeny. */
            DOSUB_SUBSET(r, sub, dmparts, ind, count, NULL);
            for (int j = 0; j < 8; j++)
               if (ci->progeny[j] != sub && ci->progeny[j] != NULL)
                   DOSUB_SUBSET(r, sub, dmparts, ind, count, ci->progeny[j]);
            
        }
        
        /* Otherwise, compute self-interaction. */
        else
           DOSELF_SUBSET(r, ci, dmparts, ind, count);
    } /* self-interaction. */
    
    /* Otherwise, it's a pair interaction. */
    else {
        
        /* Recurse? */
        if (cell_can_recurse_in_pair_dark_matter_task(ci) &&
            cell_can_recurse_in_pair_dark_matter_task(cj)) {
            
            /* Get the type of pair and flip ci/cj if needed. */
            double shift[3] = {0.0, 0.0, 0.0};
            const int sid = space_getsid(s, &ci, &cj, shift);
            
            struct cell_split_pair *csp = &cell_split_pairs[sid];
            for (int k = 0; k < csp->count; k++) {
                const int pid = csp->pairs[k].pid;
                const int pjd = csp->pairs[k].pjd;
                if (ci->progeny[pid] == sub && cj->progeny[pjd] != NULL)
                   DOSUB_SUBSET(r, ci->progeny[pid], dmparts, ind, count, cj->progeny[pjd]);
                if (ci->progeny[pid] != NULL && cj->progeny[pjd] == sub)
                   DOSUB_SUBSET(r, cj->progeny[pjd], dmparts, ind, count, ci->progeny[pid]);
            }
        }
        
        /* Otherwise, compute the pair directly. */
        else if (cell_is_active_dark_matter(ci, e) || cell_is_active_dark_matter(cj, e)) {
            
            DOPAIR_SUBSET_BRANCH(r, ci, dmparts, ind, cj, count);
        }
        
    } /* otherwise, pair interaction. */
    
    TIMER_TOC(timer_dosub_subset);
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
void DOSELF_SUBSET(struct runner *r, struct cell *restrict ci,
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
        const float hig2 = hi * hi * dm_kernel_gamma2;

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
                
                runner_iact_nonsym_dark_matter_density(r2, dx, hi, hj, pi, pj, a, H);

            }
        } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
    
    TIMER_TOC(timer_doself_subset);
}

/**
 * @brief Compute the cell self-interaction (non-symmetric).
 *
 * @param r The #runner.
 * @param c The #cell.
 */
void DOSELF1(struct runner *r, struct cell *restrict c) {
    
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    
    TIMER_TIC;
    
    struct dmpart *restrict dmparts = c->dark_matter.parts;
    const int count = c->dark_matter.count;
    
    /* Set up indt. */
    int *indt = NULL;
    int countdt = 0, firstdt = 0;
    if (posix_memalign((void **)&indt, VEC_SIZE * sizeof(int), count * sizeof(int)) != 0)
      error("Failed to allocate indt.");
    
    for (int k = 0; k < count; k++)
      if (dmpart_is_active(&dmparts[k], e)) {
        indt[countdt] = k;
        countdt += 1;
      }
    
    /* Cosmological terms */
    const float a = cosmo->a;
    const float H = cosmo->H;
    
    /* Loop over the particles in the cell. */
    for (int pid = 0; pid < count; pid++) {
        
        /* Get a pointer to the ith particle. */
        struct dmpart *restrict pi = &dmparts[pid];
        
        /* Skip inhibited particles. */
        if (dmpart_is_inhibited(pi, e)) continue;
        
        /* Get the particle position and radius. */
        double pix[3];
        for (int k = 0; k < 3; k++) pix[k] = pi->x[k];
        const float hi = pi->h;
        const float hig2 = hi * hi * dm_kernel_gamma2;
        
        /* Is the ith particle inactive? */
        if (!dmpart_is_active(pi, e)) {
            
            /* Loop over the other particles .*/
            for (int pjd = firstdt; pjd < countdt; pjd++) {
                
                /* Get a pointer to the jth particle. */
                struct dmpart *restrict pj = &dmparts[indt[pjd]];
                const float hj = pj->h;
                
#ifdef SWIFT_DEBUG_CHECKS
                /* Check that particles have been drifted to the current time */
                if (pj->ti_drift != e->ti_current)
                    error("Particle pj not drifted to current time");
#endif
                
                /* Compute the pairwise distance. */
                float r2 = 0.0f;
                float dx[3];
                for (int k = 0; k < 3; k++) {
                    dx[k] = pj->x[k] - pix[k];
                    r2 += dx[k] * dx[k];
                }
                
                /* Hit or miss? */
                if (r2 < hj * hj * dm_kernel_gamma2) {
                    
                    runner_iact_nonsym_dark_matter_density(r2, dx, hj, hi, pj, pi, a, H);
                }
            } /* loop over all other particles. */
        }
        
        /* Otherwise, interact with all candidates. */
        else {
            
            /* We caught a live one! */
            firstdt += 1;
            
            /* Loop over the other particles .*/
            for (int pjd = pid + 1; pjd < count; pjd++) {
                
                /* Get a pointer to the jth particle. */
                struct dmpart *restrict pj = &dmparts[pjd];
                
                /* Skip inhibited particles. */
                if (dmpart_is_inhibited(pj, e)) continue;
                
                const float hj = pj->h;
                
                /* Compute the pairwise distance. */
                float r2 = 0.0f;
                float dx[3];
                for (int k = 0; k < 3; k++) {
                    dx[k] = pix[k] - pj->x[k];
                    r2 += dx[k] * dx[k];
                }
                const int doj = (dmpart_is_active(pj, e)) && (r2 < hj * hj * dm_kernel_gamma2);
                const int doi = (r2 < hig2);
                
                /* Hit or miss? */
                if (doi || doj) {
                    
                    /* Which parts need to be updated? */
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
                }
            } /* loop over all other particles. */
        }
    } /* loop over all particles. */
    
    free(indt);
    
    TIMER_TOC(TIMER_DOSELF);
}

/**
 * @brief Compute the interactions within a cell (non-symmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param c The #cell.
 */
void DOSELF1_NAIVE(struct runner *r, struct cell *restrict c) {
    
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
        const float hig2 = hi * hi * dm_kernel_gamma2;
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
            const float hjg2 = hj * hj * dm_kernel_gamma2;
            const int pj_active = dmpart_is_active(pj, e);
            
#ifdef SWIFT_DEBUG_CHECKS
            /* Check that particles have been drifted to the current time */
            if (pj->ti_drift != e->ti_current)
                error("DM Particle pj not drifted to current time");
#endif
            
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
    
    TIMER_TOC(TIMER_DOSELF);
}

/**
 * @brief Determine which version of DOSELF1 needs to be called depending on the
 * optimisation level.
 *
 * @param r #runner
 * @param c #cell c
 *
 */
void DOSELF1_BRANCH(struct runner *r, struct cell *c) {
    
    const struct engine *restrict e = r->e;
    
    /* Anything to do here? */
    if (c->dark_matter.count == 0) return;
    
    /* Anything to do here? */
    if (!cell_is_active_dark_matter(c, e)) return;
    
#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
    DOSELF1_NAIVE(r, c);
#else
    DOSELF1(r, c);
#endif
}


/**
 * @brief Compute the interactions between a cell pair (non-symmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
void DOPAIR1(struct runner *r, struct cell *restrict ci, struct cell *restrict cj) {
    
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
        const float hig2 = hi * hi * dm_kernel_gamma2;
        const float pix[3] = {(float)(pi->x[0] - (cj->loc[0] + shift[0])),
            (float)(pi->x[1] - (cj->loc[1] + shift[1])),
            (float)(pi->x[2] - (cj->loc[2] + shift[2]))};
        
        /* Loop over the parts in cj. */
        for (int pjd = 0; pjd < count_j; pjd++) {
            
            /* Get a pointer to the jth particle. */
            struct dmpart *restrict pj = &dmparts_j[pjd];
            
            /* Skip inhibited particles. */
            if (dmpart_is_inhibited(pj, e)) continue;
            
            const float hj = pj->h;
            const float hjg2 = hj * hj * dm_kernel_gamma2;
            const int pj_active = dmpart_is_active(pj, e);
            
            /* Compute the pairwise distance. */
            const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                (float)(pj->x[1] - cj->loc[1]),
                (float)(pj->x[2] - cj->loc[2])};
            float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
            const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            
            /* Hit or miss? */
            if (r2 < hig2 && pi_active) {
                
                runner_iact_nonsym_dark_matter_density(r2, dx, hi, hj, pi, pj, a, H);
                
            }
            if (r2 < hjg2 && pj_active) {
                
                dx[0] = -dx[0];
                dx[1] = -dx[1];
                dx[2] = -dx[2];
                
                runner_iact_nonsym_dark_matter_density(r2, dx, hj, hi, pj, pi, a, H);
                
            }
        } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
    
    TIMER_TOC(TIMER_DOPAIR);
}


/**
 * @brief Compute the interactions between a cell pair (non-symmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
void DOPAIR1_NAIVE(struct runner *r, struct cell *restrict ci,
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
        const float hig2 = hi * hi * dm_kernel_gamma2;
        const float pix[3] = {(float)(pi->x[0] - (cj->loc[0] + shift[0])),
            (float)(pi->x[1] - (cj->loc[1] + shift[1])),
            (float)(pi->x[2] - (cj->loc[2] + shift[2]))};
        
        /* Loop over the parts in cj. */
        for (int pjd = 0; pjd < count_j; pjd++) {
            
            /* Get a pointer to the jth particle. */
            struct dmpart *restrict pj = &dmparts_j[pjd];
            
            /* Skip inhibited particles. */
            if (dmpart_is_inhibited(pj, e)) continue;
            
            const float hj = pj->h;
            const float hjg2 = hj * hj * dm_kernel_gamma2;
            const int pj_active = dmpart_is_active(pj, e);
            
            /* Compute the pairwise distance. */
            const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                (float)(pj->x[1] - cj->loc[1]),
                (float)(pj->x[2] - cj->loc[2])};
            float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
            const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            
            /* Hit or miss? */
            if (r2 < hig2 && pi_active) {
                
                runner_iact_nonsym_dark_matter_density(r2, dx, hi, hj, pi, pj, a, H);

            }
            if (r2 < hjg2 && pj_active) {
                
                dx[0] = -dx[0];
                dx[1] = -dx[1];
                dx[2] = -dx[2];
                
                runner_iact_nonsym_dark_matter_density(r2, dx, hj, hi, pj, pi, a, H);

            }
        } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
    
    TIMER_TOC(TIMER_DOPAIR);
}

/**
 * @brief Determine which version of DOPAIR1 needs to be called depending on the
 * orientation of the cells or whether DOPAIR1 needs to be called at all.
 *
 * @param r #runner
 * @param ci #cell ci
 * @param cj #cell cj
 *
 */
void DOPAIR1_BRANCH(struct runner *r, struct cell *ci, struct cell *cj) {
    
    const struct engine *restrict e = r->e;
    
    /* Anything to do here? */
    if (ci->dark_matter.count == 0 || cj->dark_matter.count == 0) return;
    
    /* Anything to do here? */
    if (!cell_is_active_dark_matter(ci, e) && !cell_is_active_dark_matter(cj, e)) return;
    
#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
    DOPAIR1_NAIVE(r, ci, cj);
#else
    DOPAIR1(r, ci, cj);
#endif
}

/**
 * @brief Compute grouped sub-cell interactions for pairs
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 *
 * @todo Hard-code the sid on the recursive calls to avoid the
 * redundant computations to find the sid on-the-fly.
 */
void DOSUB_PAIR1(struct runner *r, struct cell *ci, struct cell *cj) {
    
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
               DOSUB_PAIR1(r, ci->progeny[pid], cj->progeny[pjd]);
        }
    }
    
    /* Otherwise, compute the pair directly. */
    else if (cell_is_active_dark_matter(ci, e) || cell_is_active_dark_matter(cj, e)) {
        
        /* Compute the interactions. */
        DOPAIR1_BRANCH(r, ci, cj);
    }
    
    TIMER_TOC(TIMER_DOSUB_PAIR);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 */
void DOSUB_SELF1(struct runner *r, struct cell *ci) {
    
    TIMER_TIC;
    
    /* Should we even bother? */
    if (ci->dark_matter.count == 0 || !cell_is_active_dark_matter(ci, r->e)) return;
    
    /* Recurse? */
    if (cell_can_recurse_in_self_dark_matter_task(ci)) {
        
        /* Loop over all progeny. */
        for (int k = 0; k < 8; k++)
            if (ci->progeny[k] != NULL) {
                DOSUB_SELF1(r, ci->progeny[k]);
                for (int j = k + 1; j < 8; j++)
                   if (ci->progeny[j] != NULL)
                       DOSUB_PAIR1(r, ci->progeny[k], ci->progeny[j]);
            }
    }
    
    /* Otherwise, compute self-interaction. */
    else {
        DOSELF1_BRANCH(r, ci);
    }
    
    TIMER_TOC(TIMER_DOSUB_SELF);
}

/**
 *
 *
 *
 ********  From here SIDM specific loop functions  **********
 *
 *
 *
 **/

/**
 * @brief
 *
 * @param r The thread #runner.
 * @param ci First #cell.
 * @param cj Second #cell.
 * @param timer Are we timing this?
 */
void DOPAIR2(struct runner *r, struct cell *restrict ci, struct cell *restrict cj) {
    
    
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    const int with_cosmology = e->policy & engine_policy_cosmology;
    const struct unit_system *us = e->internal_units;
    const struct sidm_props *sidm_props = e->sidm_properties;
    
    TIMER_TIC;
    
    /* Anything to do here? */
    if (cj->dark_matter.count == 0 || ci->dark_matter.count == 0) return;
    if (!cell_is_active_dark_matter(ci, e)) return;
    
    const int count_i = ci->dark_matter.count;
    const int count_j = cj->dark_matter.count;
    struct dmpart *restrict dmparts_i = ci->dark_matter.parts;
    struct dmpart *restrict dmparts_j = cj->dark_matter.parts;
    struct sidm_history *sidm_history = &ci->dark_matter.sh;
    
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
        const integertime_t ti_step = get_integer_timestep(pi->gpart->time_bin);
        const integertime_t ti_begin = get_integer_time_begin(e->ti_current - 1, pi->gpart->time_bin);
        double dti;
        if (with_cosmology) {
            dti = cosmology_get_delta_time(e->cosmology, ti_begin,
                                           ti_begin + ti_step);
        } else {
            dti = get_timestep(pi->time_bin, e->time_base);
        }
        
        const int pi_active = dmpart_is_active(pi, e);
        const float hi = pi->h;
        const float hig2 = hi * hi * dm_kernel_gamma2;
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
                dtj = cosmology_get_delta_time(e->cosmology, ti_begin_j, ti_begin_j + ti_step_j);
            } else {
                dtj = get_timestep(pj->time_bin, e->time_base);
            }
            
            const int pj_active = dmpart_is_active(pj, e);
            const float hj = pj->h;
            const float hjg2 = hj * hj * dm_kernel_gamma2;
            
            /* Compute the pairwise distance. */
            const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                (float)(pj->x[1] - cj->loc[1]),
                (float)(pj->x[2] - cj->loc[2])};
            float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
            const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
            
            const int doi = pi_active && ((r2 < hig2) || (r2 < hjg2));
            const int doj = pj_active && ((r2 < hig2) || (r2 < hjg2));
            
            /* Hit or miss? */
            if (doi && doj) {
                
                runner_iact_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, dti, dtj, ti_begin, sidm_props, us, sidm_history, cosmo);
                
                /*runner_iact_dm_timebin(r2, dx, hi, hj, pi, pj, a, H);*/
                
            } else if (doi) {
                
                runner_iact_nonsym_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, dti, dtj, ti_begin, sidm_props, us, sidm_history, cosmo);
                
                /*runner_iact_nonsym_dm_timebin(r2, dx, hi, hj, pi, pj, a, H);*/
                
            } else if (doj) {
                
                dx[0] = -dx[0];
                dx[1] = -dx[1];
                dx[2] = -dx[2];
                
                runner_iact_nonsym_dark_matter_sidm(r2, dx, hj, hi, pj, pi, a, H, dtj, dti, ti_begin, sidm_props, us, sidm_history, cosmo);
                
                /*runner_iact_nonsym_dm_timebin(r2, dx, hj, hi, pj, pi, a, H);*/
                
            }
        } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
    
    TIMER_TOC(TIMER_DOPAIR);
}

/**
 * @brief Compute the cell self-interaction (symmetric).
 *
 * @param r The #runner.
 * @param c The #cell.
 */
void DOSELF2(struct runner *r, struct cell *restrict c) {
    
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    const int with_cosmology = e->policy & engine_policy_cosmology;
    const struct unit_system *us = e->internal_units;
    const struct sidm_props *sidm_props = e->sidm_properties;
    
    
    TIMER_TIC;
    
    struct dmpart *restrict dmparts = c->dark_matter.parts;
    const int count = c->dark_matter.count;
    struct sidm_history *sidm_history = &c->dark_matter.sh;
    
    
    /* Set up indt. */
    int *indt = NULL;
    int countdt = 0, firstdt = 0;
    if (posix_memalign((void **)&indt, VEC_SIZE * sizeof(int),
                       count * sizeof(int)) != 0)
        error("Failed to allocate indt.");
    for (int k = 0; k < count; k++)
        if (dmpart_is_active(&dmparts[k], e)) {
            indt[countdt] = k;
            countdt += 1;
        }
    
    /* Cosmological terms */
    const float a = cosmo->a;
    const float H = cosmo->H;
    
    /* Loop over the particles in the cell. */
    for (int pid = 0; pid < count; pid++) {
        
        /* Get a pointer to the ith particle. */
        struct dmpart *restrict pi = &dmparts[pid];
        
        /* Skip inhibited particles. */
        if (dmpart_is_inhibited(pi, e)) continue;
        
        /* Get i particle time-step */
        const integertime_t ti_step = get_integer_timestep(pi->gpart->time_bin);
        const integertime_t ti_begin = get_integer_time_begin(e->ti_current - 1, pi->gpart->time_bin);
        double dti;
        if (with_cosmology) {
            dti = cosmology_get_delta_time(e->cosmology, ti_begin,
                                           ti_begin + ti_step);
        } else {
            dti = get_timestep(pi->time_bin, e->time_base);
        }
        
        /* Get the particle position and radius. */
        double pix[3];
        for (int k = 0; k < 3; k++) pix[k] = pi->x[k];
        const float hi = pi->h;
        const float hig2 = hi * hi * dm_kernel_gamma2;
        
        /* Is the ith particle not active? */
        if (!dmpart_is_active(pi, e)) {
            
            /* Loop over the other particles */
            for (int pjd = firstdt; pjd < countdt; pjd++) {
                
                /* Get a pointer to the jth particle. */
                struct dmpart *restrict pj = &dmparts[indt[pjd]];
                const float hj = pj->h;
                
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
                
                /* Compute the pairwise distance. */
                float r2 = 0.0f;
                float dx[3];
                for (int k = 0; k < 3; k++) {
                    dx[k] = pj->x[k] - pix[k];
                    r2 += dx[k] * dx[k];
                }
                
                /* Hit or miss? */
                if (r2 < hig2 || r2 < hj * hj * dm_kernel_gamma2) {
                    
                    runner_iact_nonsym_dark_matter_sidm(r2, dx, hj, hi, pj, pi, a, H, dtj, dti, ti_begin, sidm_props, us, sidm_history, cosmo);
                    
                    /*runner_iact_nonsym_dm_timebin(r2, dx, hj, hi, pj, pi, a, H);*/
                }
            } /* loop over all other particles. */
        }
        
        /* Otherwise, interact with all candidates. */
        else {
            
            /* We caught a live one! */
            firstdt += 1;
            
            /* Loop over the other particles .*/
            for (int pjd = pid + 1; pjd < count; pjd++) {
                
                /* Get a pointer to the jth particle. */
                struct dmpart *restrict pj = &dmparts[pjd];
                
                /* Skip inhibited particles. */
                if (dmpart_is_inhibited(pj, e)) continue;
                
                const float hj = pj->h;
                
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
                
                /* Compute the pairwise distance. */
                float r2 = 0.0f;
                float dx[3];
                for (int k = 0; k < 3; k++) {
                    dx[k] = pix[k] - pj->x[k];
                    r2 += dx[k] * dx[k];
                }
                
                /* Hit or miss? */
                if (r2 < hig2 || r2 < hj * hj * dm_kernel_gamma2) {
                    
                    /* Does pj need to be updated too? */
                    if (dmpart_is_active(pj, e)) {
                        
                        runner_iact_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, dti, dtj, ti_begin, sidm_props, us, sidm_history, cosmo);
                        
                        /*runner_iact_dm_timebin(r2, dx, hi, hj, pi, pj, a, H);*/
                        
                    } else {
                        
                        runner_iact_nonsym_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, dti, dtj, ti_begin, sidm_props, us, sidm_history, cosmo);
                        
                        /*runner_iact_nonsym_dm_timebin(r2, dx, hi, hj, pi, pj, a, H);*/
                    }
                }
            } /* loop over all other particles. */
        }
    } /* loop over all particles. */
    
    free(indt);
    
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
void DOSELF2_NAIVE(struct runner *r, struct cell *restrict c) {
    
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    const int with_cosmology = e->policy & engine_policy_cosmology;
    const struct unit_system *us = e->internal_units;
    const struct sidm_props *sidm_props = e->sidm_properties;
    
    TIMER_TIC;
    
    /* Cosmological terms */
    const float a = cosmo->a;
    const float H = cosmo->H;
    
    const int count = c->dark_matter.count;
    struct dmpart *restrict dmparts = c->dark_matter.parts;
    struct sidm_history *sidm_history = &c->dark_matter.sh;
    
    /* Loop over the dmparts in ci. */
    for (int pid = 0; pid < count; pid++) {
        
        /* Get a hold of the ith part in ci. */
        struct dmpart *restrict pi = &dmparts[pid];
        
        /* Get i particle time-step */
        const integertime_t ti_step = get_integer_timestep(pi->gpart->time_bin);
        const integertime_t ti_begin = get_integer_time_begin(e->ti_current - 1, pi->gpart->time_bin);
        double dti;
        if (with_cosmology) {
            dti = cosmology_get_delta_time(e->cosmology, ti_begin,
                                           ti_begin + ti_step);
        } else {
            dti = get_timestep(pi->time_bin, e->time_base);
        }
        
        const int pi_active = dmpart_is_active(pi, e);
        const float hi = pi->h;
        const float hig2 = hi * hi * dm_kernel_gamma2;
        const float pix[3] = {(float)(pi->x[0] - c->loc[0]),
            (float)(pi->x[1] - c->loc[1]),
            (float)(pi->x[2] - c->loc[2])};

        /* Loop over the dmparts in cj. */
        for (int pjd = pid + 1; pjd < count; pjd++) {
            
            /* Get a pointer to the jth particle. */
            struct dmpart *restrict pj = &dmparts[pjd];
            
            const float hj = pj->h;
            const float hjg2 = hj * hj * dm_kernel_gamma2;
            const int pj_active = dmpart_is_active(pj, e);

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
                
                runner_iact_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, dti, dtj, ti_begin, sidm_props, us, sidm_history, cosmo);
                
                /*runner_iact_dm_timebin(r2, dx, hi, hj, pi, pj, a, H);*/
                
            } else if (doi) {
                
                runner_iact_nonsym_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, dti, dtj, ti_begin, sidm_props, us, sidm_history, cosmo);
                
                /*runner_iact_nonsym_dm_timebin(r2, dx, hi, hj, pi, pj, a, H);*/
                
            } else if (doj) {
                
                dx[0] = -dx[0];
                dx[1] = -dx[1];
                dx[2] = -dx[2];

                runner_iact_nonsym_dark_matter_sidm(r2, dx, hj, hi, pj, pi, a, H, dtj, dti, ti_begin, sidm_props, us, sidm_history, cosmo);
                
                /*runner_iact_nonsym_dm_timebin(r2, dx, hj, hi, pj, pi, a, H);*/
            }
        } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
    
    TIMER_TOC(TIMER_DOSELF);
}


/**
 * @brief Determine which version of DOSELF2 needs to be called depending on the
 * optimisation level.
 *
 * @param r #runner
 * @param c #cell c
 *
 */
void runner_doself2_branch_dark_matter_sidm(struct runner *r, struct cell *c) {
    
    const struct engine *restrict e = r->e;
    
    /* Anything to do here? */
    if (c->dark_matter.count == 0) return;
    
    /* Anything to do here? */
    if (!cell_is_active_dark_matter(c, e)) return;
    
    DOSELF2_NAIVE(r, c);
/*#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
    DOSELF2_NAIVE(r, c);
#else
    DOSELF2(r, c);
#endif*/
}

/**
 * @brief Compute the interactions between a cell pair (nonsymmetric case).
 *
 * Inefficient version using a brute-force algorithm.
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 */
void DOPAIR2_NAIVE(struct runner *r, struct cell *restrict ci,
                           struct cell *restrict cj) {
    
#ifdef SWIFT_DEBUG_CHECKS
    if (cj->nodeID != engine_rank) error("Should be run on a different node");
#endif
    
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    const int with_cosmology = e->policy & engine_policy_cosmology;
    const struct unit_system *us = e->internal_units;
    const struct sidm_props *sidm_props = e->sidm_properties;
    
    TIMER_TIC;
    
    const int count_i = ci->dark_matter.count;
    const int count_j = cj->dark_matter.count;
    struct dmpart *restrict dmparts_i = ci->dark_matter.parts;
    struct dmpart *restrict dmparts_j = cj->dark_matter.parts;
    struct sidm_history *sidm_history = &ci->dark_matter.sh;
    
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
        
        /* Get i particle time-step */
        const integertime_t ti_step = get_integer_timestep(pi->gpart->time_bin);
        const integertime_t ti_begin = get_integer_time_begin(e->ti_current - 1, pi->gpart->time_bin);
        double dti;
        if (with_cosmology) {
            dti = cosmology_get_delta_time(e->cosmology, ti_begin,
                                           ti_begin + ti_step);
        } else {
            dti = get_timestep(pi->time_bin, e->time_base);
        }
        
        const float hi = pi->h;
        const float hig2 = hi * hi * dm_kernel_gamma2;

        const int pi_active = dmpart_is_active(pi, e);
        const float pix[3] = {(float)(pi->x[0] - (cj->loc[0] + shift[0])),
            (float)(pi->x[1] - (cj->loc[1] + shift[1])),
            (float)(pi->x[2] - (cj->loc[2] + shift[2]))};
        
        /* Loop over the parts in cj. */
        for (int pjd = 0; pjd < count_j; pjd++) {
            
            /* Get a pointer to the jth particle. */
            struct dmpart *restrict pj = &dmparts_j[pjd];
            
            const int pj_active = dmpart_is_active(pj, e);
            const float hj = pj->h;
            const float hjg2 = hj * hj * dm_kernel_gamma2;
            
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
            
            /* Compute the pairwise distance. */
            const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                (float)(pj->x[1] - cj->loc[1]),
                (float)(pj->x[2] - cj->loc[2])};
            float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
            const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

            const int doi = pi_active && (r2 < hig2);
            const int doj = pj_active && (r2 < hjg2);
            
            /* Hit or miss? */
            if (doi && doj) {

                runner_iact_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, dti, dtj, ti_begin, sidm_props, us, sidm_history, cosmo);
                
                /*runner_iact_dm_timebin(r2, dx, hi, hj, pi, pj, a, H);*/

            } else if (doi) {
                
                runner_iact_nonsym_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, dti, dtj, ti_begin, sidm_props, us, sidm_history, cosmo);
                
                /*runner_iact_nonsym_dm_timebin(r2, dx, hi, hj, pi, pj, a, H);*/
                
            } else if (doj) {
                
                dx[0] = -dx[0];
                dx[1] = -dx[1];
                dx[2] = -dx[2];

                runner_iact_nonsym_dark_matter_sidm(r2, dx, hj, hi, pj, pi, a, H, dtj, dti, ti_begin, sidm_props, us, sidm_history, cosmo);
                
                /*runner_iact_nonsym_dm_timebin(r2, dx, hj, hi, pj, pi, a, H);*/
            }
        } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
    
    TIMER_TOC(TIMER_DOPAIR);
}


/**
 * @brief Determine which version of DOPAIR2 needs to be called depending on the
 * orientation of the cells or whether DOPAIR2 needs to be called at all.
 *
 * @param r #runner
 * @param ci #cell ci
 * @param cj #cell cj
 *
 */
void runner_dopair2_branch_dark_matter_sidm(struct runner *r, struct cell *ci, struct cell *cj) {
    
    const struct engine *restrict e = r->e;
    
    const int ci_active = cell_is_active_dark_matter(ci, e);
    const int cj_active = cell_is_active_dark_matter(cj, e);
    
    const int do_ci = (ci->dark_matter.count != 0 && ci_active);
    const int do_cj = (cj->dark_matter.count != 0 && cj_active);

    /* Anything to do here? */
    if (!do_ci && !do_cj) return;

    DOPAIR2_NAIVE(r, ci, cj);
/*#ifdef SWIFT_USE_NAIVE_INTERACTIONS
    DOPAIR2_NAIVE(r, ci, cj);
#else
    DOPAIR2(r, ci, cj);
#endif*/
}


/**
 * @brief Compute grouped sub-cell interactions for pairs (symmetric case)
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 *
 * @todo Hard-code the sid on the recursive calls to avoid the
 * redundant computations to find the sid on-the-fly.
 */
void runner_dosub_pair2_dark_matter_sidm(struct runner *r, struct cell *ci, struct cell *cj) {
    
    const struct engine *e = r->e;
    struct space *s = e->s;
    
    TIMER_TIC;
    
    /* Should we even bother? */
    const int should_do_ci = ci->dark_matter.count != 0 && cell_is_active_dark_matter(ci, e);
    const int should_do_cj = cj->dark_matter.count != 0 && cell_is_active_dark_matter(cj, e);
    
    if (!should_do_ci && !should_do_cj) return;
    
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
                runner_dosub_pair2_dark_matter_sidm(r, ci->progeny[pid], cj->progeny[pjd]);
        }
    }
    
    
    /* Otherwise, compute the pair directly. */
    else {
        /* Compute the interactions. */
        runner_dopair2_branch_dark_matter_sidm(r, ci, cj);
    }
    
    TIMER_TOC(TIMER_DOSUB_PAIR);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks (symmetric case)
 *
 * @param r The #runner.
 * @param ci The first #cell.
 */
void runner_dosub_self2_dark_matter_sidm(struct runner *r, struct cell *ci) {
    
    TIMER_TIC;
    
    /* Should we even bother? */
    if (ci->dark_matter.count == 0 || !cell_is_active_dark_matter(ci, r->e)) return;
    
    /* Recurse? */
    if (cell_can_recurse_in_self_dark_matter_task(ci)) {
        
        /* Loop over all progeny. */
        for (int k = 0; k < 8; k++)
            if (ci->progeny[k] != NULL) {
                runner_dosub_self2_dark_matter_sidm(r, ci->progeny[k]);
                for (int j = k + 1; j < 8; j++)
                    if (ci->progeny[j] != NULL)
                        runner_dosub_pair2_dark_matter_sidm(r, ci->progeny[k], ci->progeny[j]);
            }
    }
    
    /* Otherwise, compute self-interaction. */
    else {
        runner_doself2_branch_dark_matter_sidm(r, ci);
    }
    
    TIMER_TOC(TIMER_DOSUB_SELF);
}



