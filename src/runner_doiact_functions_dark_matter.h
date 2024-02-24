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

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
                /* Check that particles have been drifted to the current time */
                if (pi->ti_drift != e->ti_current)
                    error("Particle pi not drifted to current time");
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
                
#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
                /* Check that particles have been drifted to the current time */
                if (pi->ti_drift != e->ti_current)
                    error("Particle pi not drifted to current time");
                if (pj->ti_drift != e->ti_current)
                    error("Particle pj not drifted to current time");
#endif
                
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
            
#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
            /* Check that particles have been drifted to the current time */
            if (pi->ti_drift != e->ti_current)
                error("Particle pi not drifted to current time");
            if (pj->ti_drift != e->ti_current)
                error("Particle pj not drifted to current time");
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

    /* Did we mess up the recursion? */
    if (c->dark_matter.h_max_old * kernel_gamma > c->dmin)
        error("Cell smaller than smoothing length");

    /* Check that cells are drifted. */
    if (!cell_are_dmpart_drifted(c, e)) error("Interacting undrifted cell.");

#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
    DOSELF1_NAIVE(r, c);
#else
    DOSELF1(r, c);
#endif
}


/**
 * @brief Compute the interactions between a cell pair (non-symmetric case).
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param sid The direction of the pair.
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR1(struct runner *r, struct cell *restrict ci, struct cell *restrict cj,
             const int sid, const double *shift) {

    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    
    TIMER_TIC;

    /* Get the cutoff shift. */
    double rshift = 0.0;
    for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

    /* Pick-out the sorted lists. */
    const struct sort_entry *restrict sort_i = cell_get_dark_matter_sorts(ci, sid);
    const struct sort_entry *restrict sort_j = cell_get_dark_matter_sorts(cj, sid);

#ifdef SWIFT_DEBUG_CHECKS
    /* Some constants used to checks that the parts are in the right frame */
  const float shift_threshold_x = 2.02 * ci->width[0] +
      2.02 * max(ci->dark_matter.dx_max_part, cj->dark_matter.dx_max_part);
  const float shift_threshold_y = 2.02 * ci->width[1] +
      2.02 * max(ci->dark_matter.dx_max_part, cj->dark_matter.dx_max_part);
  const float shift_threshold_z = 2.02 * ci->width[2] +
      2.02 * max(ci->dark_matter.dx_max_part, cj->dark_matter.dx_max_part);
#endif /* SWIFT_DEBUG_CHECKS */

    /* Get some other useful values. */
    const double hi_max = ci->dark_matter.h_max * dm_kernel_gamma - rshift;
    const double hj_max = cj->dark_matter.h_max * dm_kernel_gamma;
    const int count_i = ci->dark_matter.count;
    const int count_j = cj->dark_matter.count;
    struct dmpart *restrict parts_i = ci->dark_matter.parts;
    struct dmpart *restrict parts_j = cj->dark_matter.parts;
    const double di_max = sort_i[count_i - 1].d - rshift;
    const double dj_min = sort_j[0].d;
    const float dx_max = (ci->dark_matter.dx_max_sort + cj->dark_matter.dx_max_sort);

    /* Cosmological terms and physical constants */
    const float a = cosmo->a;
    const float H = cosmo->H;

    if (cell_is_active_dark_matter(ci, e)) {

        /* Loop over the parts in ci. */
        for (int pid = count_i - 1;
             pid >= 0 && sort_i[pid].d + hi_max + dx_max > dj_min; pid--) {

            /* Get a hold of the ith part in ci. */
            struct dmpart *restrict pi = &parts_i[sort_i[pid].i];
            const float hi = pi->h;

            /* Skip inactive particles */
            if (!dmpart_is_active(pi, e)) continue;

            /* Is there anything we need to interact with ? */
            const double di = sort_i[pid].d + hi * dm_kernel_gamma + dx_max - rshift;
            if (di < dj_min) continue;

            /* Get some additional information about pi */
            const float hig2 = hi * hi * dm_kernel_gamma2;
            const float pix = pi->x[0] - (cj->loc[0] + shift[0]);
            const float piy = pi->x[1] - (cj->loc[1] + shift[1]);
            const float piz = pi->x[2] - (cj->loc[2] + shift[2]);

            /* Loop over the parts in cj. */
            for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

                /* Recover pj */
                struct dmpart *pj = &parts_j[sort_j[pjd].i];

                /* Skip inhibited particles. */
                if (dmpart_is_inhibited(pj, e)) continue;

                const float hj = pj->h;
                const float pjx = pj->x[0] - cj->loc[0];
                const float pjy = pj->x[1] - cj->loc[1];
                const float pjz = pj->x[2] - cj->loc[2];

                /* Compute the pairwise distance. */
                float dx[3] = {pix - pjx, piy - pjy, piz - pjz};
                const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
                /* Check that particles are in the correct frame after the shifts */
                if (pix > shift_threshold_x || pix < -shift_threshold_x)
                  error(
                      "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
                      pix, ci->width[0]);
                if (piy > shift_threshold_y || piy < -shift_threshold_y)
                  error(
                      "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
                      piy, ci->width[1]);
                if (piz > shift_threshold_z || piz < -shift_threshold_z)
                  error(
                      "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
                      piz, ci->width[2]);
                if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
                  error(
                      "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
                      pjx, ci->width[0]);
                if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
                  error(
                      "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
                      pjy, ci->width[1]);
                if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
                  error(
                      "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
                      pjz, ci->width[2]);

#if defined(DO_DRIFT_DEBUG_CHECKS)
                /* Check that particles have been drifted to the current time */
                if (pi->ti_drift != e->ti_current)
                  error("Particle pi not drifted to current time");
                if (pj->ti_drift != e->ti_current)
                  error("Particle pj not drifted to current time");
#endif
#endif

                /* Hit or miss? */
                if (r2 < hig2) {

                    runner_iact_nonsym_dark_matter_density(r2, dx, hi, hj, pi, pj, a, H);

                }
            } /* loop over the parts in cj. */
        }   /* loop over the parts in ci. */
    }     /* Cell ci is active */

    if (cell_is_active_dark_matter(cj, e)) {

        /* Loop over the parts in cj. */
        for (int pjd = 0; pjd < count_j && sort_j[pjd].d - hj_max - dx_max < di_max;
             pjd++) {

            /* Get a hold of the jth part in cj. */
            struct dmpart *pj = &parts_j[sort_j[pjd].i];
            const float hj = pj->h;

            /* Skip inactive particles */
            if (!dmpart_is_active(pj, e)) continue;

            /* Is there anything we need to interact with ? */
            const double dj = sort_j[pjd].d - hj * dm_kernel_gamma - dx_max + rshift;
            if (dj - rshift > di_max) continue;

            /* Get some additional information about pj */
            const float hjg2 = hj * hj * dm_kernel_gamma2;
            const float pjx = pj->x[0] - cj->loc[0];
            const float pjy = pj->x[1] - cj->loc[1];
            const float pjz = pj->x[2] - cj->loc[2];

            /* Loop over the parts in ci. */
            for (int pid = count_i - 1; pid >= 0 && sort_i[pid].d > dj; pid--) {

                /* Recover pi */
                struct dmpart *pi = &parts_i[sort_i[pid].i];

                /* Skip inhibited particles. */
                if (dmpart_is_inhibited(pi, e)) continue;

                const float hi = pi->h;
                const float pix = pi->x[0] - (cj->loc[0] + shift[0]);
                const float piy = pi->x[1] - (cj->loc[1] + shift[1]);
                const float piz = pi->x[2] - (cj->loc[2] + shift[2]);

                /* Compute the pairwise distance. */
                float dx[3] = {pjx - pix, pjy - piy, pjz - piz};
                const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
                /* Check that particles are in the correct frame after the shifts */
                if (pix > shift_threshold_x || pix < -shift_threshold_x)
                  error(
                      "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
                      pix, ci->width[0]);
                if (piy > shift_threshold_y || piy < -shift_threshold_y)
                  error(
                      "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
                      piy, ci->width[1]);
                if (piz > shift_threshold_z || piz < -shift_threshold_z)
                  error(
                      "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
                      piz, ci->width[2]);
                if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
                  error(
                      "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
                      pjx, ci->width[0]);
                if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
                  error(
                      "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
                      pjy, ci->width[1]);
                if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
                  error(
                      "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
                      pjz, ci->width[2]);

#if defined(DO_DRIFT_DEBUG_CHECKS)
                /* Check that particles have been drifted to the current time */
                if (pi->ti_drift != e->ti_current)
                  error("Particle pi not drifted to current time");
                if (pj->ti_drift != e->ti_current)
                  error("Particle pj not drifted to current time");
#endif
#endif

                /* Hit or miss? */
                if (r2 < hjg2) {

                    runner_iact_nonsym_dark_matter_density(r2, dx, hj, hi, pj, pi, a, H);

                }
            } /* loop over the parts in ci. */
        }     /* loop over the parts in cj. */
    }         /* Cell cj is active */
    
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
            
#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
            /* Check that particles have been drifted to the current time */
            if (pi->ti_drift != e->ti_current)
                error("Particle pi not drifted to current time");
            if (pj->ti_drift != e->ti_current)
                error("Particle pj not drifted to current time");
#endif

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
    }     /* loop over the parts in ci. */
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

    /* Check that cells are drifted. */
    if (!cell_are_dmpart_drifted(ci, e) || !cell_are_dmpart_drifted(cj, e))
        error("Interacting undrifted cells.");

    /* Get the sort ID. */
    double shift[3] = {0.0, 0.0, 0.0};
    const int sid = space_getsid(e->s, &ci, &cj, shift);

    /* Have the cells been sorted? */
    if (!(ci->dark_matter.sorted & (1 << sid)) ||
        ci->dark_matter.dx_max_sort_old > space_maxreldx * ci->dmin)
        error("Interacting unsorted cells.");
    if (!(cj->dark_matter.sorted & (1 << sid)) ||
        cj->dark_matter.dx_max_sort_old > space_maxreldx * cj->dmin)
        error("Interacting unsorted cells.");

#ifdef SWIFT_DEBUG_CHECKS
        /* Pick-out the sorted lists. */
  const struct sort_entry *restrict sort_i = cell_get_dark_matter_sorts(ci, sid);
  const struct sort_entry *restrict sort_j = cell_get_dark_matter_sorts(cj, sid);

  /* Check that the dx_max_sort values in the cell are indeed an upper
     bound on particle movement. */
  for (int pid = 0; pid < ci->dark_matter.count; pid++) {
    const struct dmpart *p = &ci->dark_matter.parts[sort_i[pid].i];
    if (dmpart_is_inhibited(p, e)) continue;

    const float d = p->x[0] * runner_shift[sid][0] +
                    p->x[1] * runner_shift[sid][1] +
                    p->x[2] * runner_shift[sid][2];
    if (fabsf(d - sort_i[pid].d) - ci->dark_matter.dx_max_sort >
            1.0e-4 * max(fabsf(d), ci->dark_matter.dx_max_sort_old) &&
        fabsf(d - sort_i[pid].d) - ci->dark_matter.dx_max_sort >
            ci->width[0] * 1.0e-10)
      error(
          "particle shift diff exceeds dx_max_sort in cell ci. ci->nodeID=%d "
          "cj->nodeID=%d d=%e sort_i[pid].d=%e ci->dark_matter.dx_max_sort=%e "
          "ci->dark_matter.dx_max_sort_old=%e",
          ci->nodeID, cj->nodeID, d, sort_i[pid].d, ci->dark_matter.dx_max_sort,
          ci->dark_matter.dx_max_sort_old);
  }
  for (int pjd = 0; pjd < cj->dark_matter.count; pjd++) {
    const struct dmpart *p = &cj->dark_matter.parts[sort_j[pjd].i];
    if (dmpart_is_inhibited(p, e)) continue;

    const float d = p->x[0] * runner_shift[sid][0] +
                    p->x[1] * runner_shift[sid][1] +
                    p->x[2] * runner_shift[sid][2];
    if ((fabsf(d - sort_j[pjd].d) - cj->dark_matter.dx_max_sort) >
            1.0e-4 * max(fabsf(d), cj->dark_matter.dx_max_sort_old) &&
        (fabsf(d - sort_j[pjd].d) - cj->dark_matter.dx_max_sort) >
            cj->width[0] * 1.0e-10)
      error(
          "particle shift diff exceeds dx_max_sort in cell cj. cj->nodeID=%d "
          "ci->nodeID=%d d=%e sort_j[pjd].d=%e cj->dark_matter.dx_max_sort=%e "
          "cj->dark_matter.dx_max_sort_old=%e",
          cj->nodeID, ci->nodeID, d, sort_j[pjd].d, cj->dark_matter.dx_max_sort,
          cj->dark_matter.dx_max_sort_old);
  }
#endif /* SWIFT_DEBUG_CHECKS */

#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
    DOPAIR1_NAIVE(r, ci, cj);
#else
    DOPAIR1(r, ci, cj, sid, shift);
#endif
}

/**
 * @brief Compute grouped sub-cell interactions for pairs
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param gettimer Do we have a timer ?
 *
 * @todo Hard-code the sid on the recursive calls to avoid the
 * redundant computations to find the sid on-the-fly.
 */
void DOSUB_PAIR1(struct runner *r, struct cell *ci, struct cell *cj,
                 int gettimer) {
    
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
               DOSUB_PAIR1(r, ci->progeny[pid], cj->progeny[pjd], 0);
        }
    }
    
    /* Otherwise, compute the pair directly. */
    else if (cell_is_active_dark_matter(ci, e) || cell_is_active_dark_matter(cj, e)) {

        /* Make sure both cells are drifted to the current timestep. */
        if (!cell_are_dmpart_drifted(ci, e) || !cell_are_dmpart_drifted(cj, e))
            error("Interacting undrifted cells.");

        /* Do any of the cells need to be sorted first? */
        if (!(ci->dark_matter.sorted & (1 << sid)) ||
            ci->dark_matter.dx_max_sort_old > ci->dmin * space_maxreldx)
            error(
                "Interacting unsorted cell. ci->dark_matter.dx_max_sort_old=%e ci->dmin=%e "
                "ci->sorted=%d sid=%d",
                ci->dark_matter.dx_max_sort_old, ci->dmin, ci->dark_matter.sorted, sid);
        if (!(cj->dark_matter.sorted & (1 << sid)) ||
            cj->dark_matter.dx_max_sort_old > cj->dmin * space_maxreldx)
            error(
                "Interacting unsorted cell. cj->dark_matter.dx_max_sort_old=%e cj->dmin=%e "
                "cj->sorted=%d sid=%d",
                cj->dark_matter.dx_max_sort_old, cj->dmin, cj->dark_matter.sorted, sid);

        /* Compute the interactions. */
        DOPAIR1_BRANCH(r, ci, cj);
    }
    if (gettimer) TIMER_TOC(TIMER_DOSUB_PAIR);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF1(struct runner *r, struct cell *ci, int gettimer) {
    
    TIMER_TIC;
    
    /* Should we even bother? */
    if (ci->dark_matter.count == 0 || !cell_is_active_dark_matter(ci, r->e)) return;
    
    /* Recurse? */
    if (cell_can_recurse_in_self_dark_matter_task(ci)) {
        
        /* Loop over all progeny. */
        for (int k = 0; k < 8; k++)
            if (ci->progeny[k] != NULL) {
                DOSUB_SELF1(r, ci->progeny[k], 0);
                for (int j = k + 1; j < 8; j++)
                   if (ci->progeny[j] != NULL)
                       DOSUB_PAIR1(r, ci->progeny[k], ci->progeny[j], 0);
            }
    }
    
    /* Otherwise, compute self-interaction. */
    else {

        /* Drift the cell to the current timestep if needed. */
        if (!cell_are_dmpart_drifted(ci, r->e)) error("Interacting undrifted cell.");

        DOSELF1_BRANCH(r, ci);
    }
    if (gettimer) TIMER_TOC(TIMER_DOSUB_SELF);
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
 * @brief Compute the interactions between a cell pair (symmetric)
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param sid The direction of the pair
 * @param shift The shift vector to apply to the particles in ci.
 */
void DOPAIR2(struct runner *r, struct cell *ci, struct cell *cj, const int sid,
             const double *shift) {
    
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    const struct unit_system *us = e->internal_units;
    const struct sidm_props *sidm_props = e->sidm_properties;

    const int with_cosmology = e->policy & engine_policy_cosmology;
    const double time_base = e->time_base;
    const integertime_t t_current = e->ti_current;

    TIMER_TIC;

    /* Get the cutoff shift. */
    double rshift = 0.0;
    for (int k = 0; k < 3; k++) rshift += shift[k] * runner_shift[sid][k];

    /* Pick-out the sorted lists. */
    struct sort_entry *restrict sort_i = cell_get_dark_matter_sorts(ci, sid);
    struct sort_entry *restrict sort_j = cell_get_dark_matter_sorts(cj, sid);

#ifdef SWIFT_DEBUG_CHECKS
    /* Some constants used to checks that the parts are in the right frame */
  const float shift_threshold_x = 2.02 * ci->width[0] +
      2.02 * max(ci->dark_matter.dx_max_part, cj->dark_matter.dx_max_part);
  const float shift_threshold_y = 2.02 * ci->width[1] +
      2.02 * max(ci->dark_matter.dx_max_part, cj->dark_matter.dx_max_part);
  const float shift_threshold_z = 2.02 * ci->width[2] +
      2.02 * max(ci->dark_matter.dx_max_part, cj->dark_matter.dx_max_part);
#endif /* SWIFT_DEBUG_CHECKS */

    /* Get some other useful values. */
    const double hi_max = ci->dark_matter.h_max;
    const double hj_max = cj->dark_matter.h_max;
    const int count_i = ci->dark_matter.count;
    const int count_j = cj->dark_matter.count;
    struct dmpart *restrict parts_i = ci->dark_matter.parts;
    struct dmpart *restrict parts_j = cj->dark_matter.parts;

    /* Cosmological terms and physical constants */
    const float a = cosmo->a;
    const float H = cosmo->H;

    /* Maximal displacement since last rebuild */
    const double dx_max = (ci->dark_matter.dx_max_sort + cj->dark_matter.dx_max_sort);

    /* Position on the axis of the particles closest to the interface */
    const double di_max = sort_i[count_i - 1].d;
    const double dj_min = sort_j[0].d;

    /* Shifts to apply to the particles to be in a good frame */
    const double shift_i[3] = {cj->loc[0] + shift[0], cj->loc[1] + shift[1],
                               cj->loc[2] + shift[2]};
    const double shift_j[3] = {cj->loc[0], cj->loc[1], cj->loc[2]};

    int count_active_i = 0, count_active_j = 0;
    struct sort_entry *restrict sort_active_i = NULL;
    struct sort_entry *restrict sort_active_j = NULL;

    // MATTHIEU: temporary disable this optimization
    if (0 /*&& cell_is_all_active_dark_matter(cj, e)*/) {
        /* If everybody is active don't bother copying */
        sort_active_i = sort_i;
        count_active_i = count_i;
    } else if (cell_is_active_dark_matter(ci, e)) {
        if (posix_memalign((void **)&sort_active_i, SWIFT_CACHE_ALIGNMENT,
                           sizeof(struct sort_entry) * count_i) != 0)
            error("Failed to allocate active sortlists.");

        /* Collect the active particles in ci */
        for (int k = 0; k < count_i; k++) {
            if (dmpart_is_active(&parts_i[sort_i[k].i], e)) {
                sort_active_i[count_active_i] = sort_i[k];
                count_active_i++;
            }
        }
    }

    // MATTHIEU: temporary disable this optimization
    if (0 /*&& cell_is_all_active_dark_matter(cj, e)*/) {
        /* If everybody is active don't bother copying */
        sort_active_j = sort_j;
        count_active_j = count_j;
    } else if (cell_is_active_dark_matter(cj, e)) {
        if (posix_memalign((void **)&sort_active_j, SWIFT_CACHE_ALIGNMENT,
                           sizeof(struct sort_entry) * count_j) != 0)
            error("Failed to allocate active sortlists.");

        /* Collect the active dmparticles in cj */
        for (int k = 0; k < count_j; k++) {
            if (dmpart_is_active(&parts_j[sort_j[k].i], e)) {
                sort_active_j[count_active_j] = sort_j[k];
                count_active_j++;
            }
        }
    }

    /* Loop over *all* the parts in ci starting from the centre until
       we are out of range of anything in cj (using the maximal hi). */
    for (int pid = count_i - 1;
         pid >= 0 &&
         sort_i[pid].d + hi_max * dm_kernel_gamma + dx_max - rshift > dj_min;
         pid--) {

        /* Get a hold of the ith part in ci. */
        struct dmpart *pi = &parts_i[sort_i[pid].i];

        /* Skip inhibited particles. */
        if (dmpart_is_inhibited(pi, e)) continue;

        const float hi = pi->h;

        /* Is there anything we need to interact with (for this specific hi) ? */
        const double di = sort_i[pid].d + hi * dm_kernel_gamma + dx_max - rshift;
        if (di < dj_min) continue;

        /* Get some additional information about pi */
        const float hig2 = hi * hi * dm_kernel_gamma2;
        const float pix = pi->x[0] - shift_i[0];
        const float piy = pi->x[1] - shift_i[1];
        const float piz = pi->x[2] - shift_i[2];

        /* Do we need to only check active parts in cj
           (i.e. pi does not need updating) ? */
        if (!dmpart_is_active(pi, e)) {

            /* Loop over the *active* parts in cj within range of pi */
            for (int pjd = 0; pjd < count_active_j && sort_active_j[pjd].d < di;
                 pjd++) {

                /* Recover pj */
                struct dmpart *pj = &parts_j[sort_active_j[pjd].i];

                /* Skip inhibited particles. */
                if (dmpart_is_inhibited(pj, e)) continue;

                const float hj = pj->h;

                /* Get the position of pj in the right frame */
                const float pjx = pj->x[0] - shift_j[0];
                const float pjy = pj->x[1] - shift_j[1];
                const float pjz = pj->x[2] - shift_j[2];

                /* Compute the pairwise distance. */
                const float dx[3] = {pjx - pix, pjy - piy, pjz - piz};
                const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
                /* Check that particles are in the correct frame after the shifts */
                if (pix > shift_threshold_x || pix < -shift_threshold_x)
                  error(
                      "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
                      pix, ci->width[0]);
                if (piy > shift_threshold_y || piy < -shift_threshold_y)
                  error(
                      "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
                      piy, ci->width[1]);
                if (piz > shift_threshold_z || piz < -shift_threshold_z)
                  error(
                      "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
                      piz, ci->width[2]);
                if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
                  error(
                      "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
                      pjx, ci->width[0]);
                if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
                  error(
                      "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
                      pjy, ci->width[1]);
                if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
                  error(
                      "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
                      pjz, ci->width[2]);

#if defined(DO_DRIFT_DEBUG_CHECKS)
                /* Check that particles have been drifted to the current time */
                if (pi->ti_drift != e->ti_current)
                  error("Particle pi not drifted to current time");

                if (pj->ti_drift != e->ti_current)
                  error("Particle pj not drifted to current time");
#endif
#endif
                /* Hit or miss?
                   (note that we will do the other condition in the reverse loop) */
                if (r2 < hig2) {

                    runner_iact_nonsym_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                                        t_current, cosmo, with_cosmology, sidm_props, us);

                }
            } /* loop over the parts in cj. */
        }   /* loop over the parts in ci. */

        else { /* pi is active, we may need to update pi and pj */

            /* Loop over *all* the parts in cj in range of pi. */
            for (int pjd = 0; pjd < count_j && sort_j[pjd].d < di; pjd++) {

                /* Recover pj */
                struct dmpart *pj = &parts_j[sort_j[pjd].i];

                /* Skip inhibited particles. */
                if (dmpart_is_inhibited(pj, e)) continue;

                const float hj = pj->h;

                /* Get the position of pj in the right frame */
                const float pjx = pj->x[0] - shift_j[0];
                const float pjy = pj->x[1] - shift_j[1];
                const float pjz = pj->x[2] - shift_j[2];

                /* Compute the pairwise distance. */
                const float dx[3] = {pix - pjx, piy - pjy, piz - pjz};
                const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
                /* Check that particles are in the correct frame after the shifts */
                if (pix > shift_threshold_x || pix < -shift_threshold_x)
                  error(
                      "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
                      pix, ci->width[0]);
                if (piy > shift_threshold_y || piy < -shift_threshold_y)
                  error(
                      "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
                      piy, ci->width[1]);
                if (piz > shift_threshold_z || piz < -shift_threshold_z)
                  error(
                      "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
                      piz, ci->width[2]);
                if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
                  error(
                      "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
                      pjx, ci->width[0]);
                if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
                  error(
                      "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
                      pjy, ci->width[1]);
                if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
                  error(
                      "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
                      pjz, ci->width[2]);

#if defined(DO_DRIFT_DEBUG_CHECKS)
                /* Check that particles have been drifted to the current time */
                if (pi->ti_drift != e->ti_current)
                  error("Particle pi not drifted to current time");

                if (pj->ti_drift != e->ti_current)
                  error("Particle pj not drifted to current time");
#endif
#endif
                /* Hit or miss?
                   (note that we will do the other condition in the reverse loop) */
                if (r2 < hig2) {

                    /* Does pj need to be updated too? */
                    if (dmpart_is_active(pj, e)) {

                        runner_iact_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                                     t_current, cosmo, with_cosmology, sidm_props, us);

                    } else {

                        runner_iact_nonsym_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                                            t_current, cosmo, with_cosmology, sidm_props, us);
                    }
                }
            } /* loop over the parts in cj. */
        }     /* Is pi active? */
    }         /* Loop over all ci */

    /* Loop over *all* the parts in cj starting from the centre until
       we are out of range of anything in ci (using the maximal hj). */
    for (int pjd = 0;
         pjd < count_j &&
         sort_j[pjd].d - hj_max * dm_kernel_gamma - dx_max < di_max - rshift;
         pjd++) {

        /* Get a hold of the jth part in cj. */
        struct dmpart *pj = &parts_j[sort_j[pjd].i];

        /* Skip inhibited particles. */
        if (dmpart_is_inhibited(pj, e)) continue;

        const float hj = pj->h;

        /* Is there anything we need to interact with (for this specific hj) ? */
        const double dj = sort_j[pjd].d - hj * dm_kernel_gamma - dx_max;
        if (dj > di_max - rshift) continue;

        /* Get some additional information about pj */
        const float hjg2 = hj * hj * dm_kernel_gamma2;
        const float pjx = pj->x[0] - shift_j[0];
        const float pjy = pj->x[1] - shift_j[1];
        const float pjz = pj->x[2] - shift_j[2];

        /* Do we need to only check active parts in ci
           (i.e. pj does not need updating) ? */
        if (!dmpart_is_active(pj, e)) {

            /* Loop over the *active* parts in ci. */
            for (int pid = count_active_i - 1;
                 pid >= 0 && sort_active_i[pid].d - rshift > dj; pid--) {

                /* Recover pi */
                struct dmpart *pi = &parts_i[sort_active_i[pid].i];

                /* Skip inhibited particles */
                if (dmpart_is_inhibited(pi, e)) continue;

                const float hi = pi->h;
                const float hig2 = hi * hi * dm_kernel_gamma2;

                /* Get the position of pi in the right frame */
                const float pix = pi->x[0] - shift_i[0];
                const float piy = pi->x[1] - shift_i[1];
                const float piz = pi->x[2] - shift_i[2];

                /* Compute the pairwise distance. */
                const float dx[3] = {pix - pjx, piy - pjy, piz - pjz};
                const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
                /* Check that particles are in the correct frame after the shifts */
                if (pix > shift_threshold_x || pix < -shift_threshold_x)
                  error(
                      "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
                      pix, ci->width[0]);
                if (piy > shift_threshold_y || piy < -shift_threshold_y)
                  error(
                      "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
                      piy, ci->width[1]);
                if (piz > shift_threshold_z || piz < -shift_threshold_z)
                  error(
                      "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
                      piz, ci->width[2]);
                if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
                  error(
                      "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
                      pjx, ci->width[0]);
                if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
                  error(
                      "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
                      pjy, ci->width[1]);
                if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
                  error(
                      "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
                      pjz, ci->width[2]);

#if defined(DO_DRIFT_DEBUG_CHECKS)
                /* Check that particles have been drifted to the current time */
                if (pi->ti_drift != e->ti_current)
                  error("Particle pi not drifted to current time");
                if (pj->ti_drift != e->ti_current)
                  error("Particle pj not drifted to current time");
#endif
#endif
                /* Hit or miss?
                   (note that we must avoid the r2 < hig2 cases we already processed) */
                if (r2 < hjg2 && r2 >= hig2) {

                    runner_iact_nonsym_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                                        t_current, cosmo, with_cosmology, sidm_props, us);
                }
            } /* loop over the active parts in ci. */
        }

        else { /* pj is active, we may need to update pj and pi */

            /* Loop over *all* the parts in ci. */
            for (int pid = count_i - 1; pid >= 0 && sort_i[pid].d - rshift > dj;
                 pid--) {

                /* Recover pi */
                struct dmpart *pi = &parts_i[sort_i[pid].i];

                /* Skip inhibited particles. */
                if (dmpart_is_inhibited(pi, e)) continue;

                const float hi = pi->h;
                const float hig2 = hi * hi * dm_kernel_gamma2;

                /* Get the position of pi in the right frame */
                const float pix = pi->x[0] - shift_i[0];
                const float piy = pi->x[1] - shift_i[1];
                const float piz = pi->x[2] - shift_i[2];

                /* Compute the pairwise distance. */
                const float dx[3] = {pjx - pix, pjy - piy, pjz - piz};
                const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#ifdef SWIFT_DEBUG_CHECKS
                        /* Check that particles are in the correct frame after the shifts */
                if (pix > shift_threshold_x || pix < -shift_threshold_x)
                  error(
                      "Invalid particle position in X for pi (pix=%e ci->width[0]=%e)",
                      pix, ci->width[0]);
                if (piy > shift_threshold_y || piy < -shift_threshold_y)
                  error(
                      "Invalid particle position in Y for pi (piy=%e ci->width[1]=%e)",
                      piy, ci->width[1]);
                if (piz > shift_threshold_z || piz < -shift_threshold_z)
                  error(
                      "Invalid particle position in Z for pi (piz=%e ci->width[2]=%e)",
                      piz, ci->width[2]);
                if (pjx > shift_threshold_x || pjx < -shift_threshold_x)
                  error(
                      "Invalid particle position in X for pj (pjx=%e ci->width[0]=%e)",
                      pjx, ci->width[0]);
                if (pjy > shift_threshold_y || pjy < -shift_threshold_y)
                  error(
                      "Invalid particle position in Y for pj (pjy=%e ci->width[1]=%e)",
                      pjy, ci->width[1]);
                if (pjz > shift_threshold_z || pjz < -shift_threshold_z)
                  error(
                      "Invalid particle position in Z for pj (pjz=%e ci->width[2]=%e)",
                      pjz, ci->width[2]);

#if defined(DO_DRIFT_DEBUG_CHECKS)
                /* Check that particles have been drifted to the current time */
                if (pi->ti_drift != e->ti_current)
                  error("Particle pi not drifted to current time");
                if (pj->ti_drift != e->ti_current)
                  error("Particle pj not drifted to current time");
#endif
#endif
                /* Hit or miss?
                   (note that we must avoid the r2 < hig2 cases we already processed) */
                if (r2 < hjg2 && r2 >= hig2) {

                    /* Does pi need to be updated too? */
                    if (dmpart_is_active(pi, e)) {

                        runner_iact_dark_matter_sidm(r2, dx, hj, hi, pj, pi, a, H, time_base,
                                                     t_current, cosmo, with_cosmology, sidm_props, us);
                    } else {

                        runner_iact_nonsym_dark_matter_sidm(r2, dx, hj, hi, pj, pi, a, H, time_base,
                                                            t_current, cosmo, with_cosmology, sidm_props, us);
                    }
                }
            } /* loop over the parts in ci. */
        }   /* Is pj active? */
    }     /* Loop over all cj */

    /* Clean-up if necessary */
    if (cell_is_active_dark_matter(ci, e))
        free(sort_active_i);
    if (cell_is_active_dark_matter(cj, e))
        free(sort_active_j);

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
void DOPAIR2_NAIVE(struct runner *r, struct cell *restrict ci,
                   struct cell *restrict cj) {

    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    const struct unit_system *us = e->internal_units;
    const struct sidm_props *sidm_props = e->sidm_properties;
    const int with_cosmology = (e->policy & engine_policy_cosmology);
    const double time_base = e->time_base;
    const integertime_t t_current = e->ti_current;

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

            const int pj_active = dmpart_is_active(pj, e);
            const float hj = pj->h;
            const float hjg2 = hj * hj * dm_kernel_gamma2;

            /* Compute the pairwise distance. */
            const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                                  (float)(pj->x[1] - cj->loc[1]),
                                  (float)(pj->x[2] - cj->loc[2])};
            float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
            const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
            /* Check that particles have been drifted to the current time */
              if (pi->ti_drift != e->ti_current)
                error("Particle pi not drifted to current time");
              if (pj->ti_drift != e->ti_current)
                error("Particle pj not drifted to current time");
#endif

            /* Hit or miss? */
            if (r2 < hig2 || r2 < hjg2) {

                if (pi_active && pj_active) {

                    runner_iact_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                                 t_current, cosmo, with_cosmology, sidm_props, us);
                } else if (pi_active) {

                    runner_iact_nonsym_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                                        t_current, cosmo, with_cosmology, sidm_props, us);
                } else if (pj_active) {

                    dx[0] = -dx[0];
                    dx[1] = -dx[1];
                    dx[2] = -dx[2];

                    runner_iact_nonsym_dark_matter_sidm(r2, dx, hj, hi, pj, pi, a, H, time_base,
                                                        t_current, cosmo, with_cosmology, sidm_props, us);
                }
            }
        } /* loop over the parts in cj. */
    }     /* loop over the parts in ci. */
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
void DOPAIR2_BRANCH(struct runner *r, struct cell *ci, struct cell *cj) {

    const struct engine *restrict e = r->e;

    /* Anything to do here? */
    if (ci->dark_matter.count == 0 || cj->dark_matter.count == 0) return;

    /* Anything to do here? */
    if (!cell_is_active_dark_matter(ci, e) && !cell_is_active_dark_matter(cj, e)) return;

    /* Check that cells are drifted. */
    if (!cell_are_dmpart_drifted(ci, e) || !cell_are_dmpart_drifted(cj, e))
        error("Interacting undrifted cells.");

    /* Get the sort ID. */
    double shift[3] = {0.0, 0.0, 0.0};
    const int sid = space_getsid(e->s, &ci, &cj, shift);

    /* Have the cells been sorted? */
    if (!(ci->dark_matter.sorted & (1 << sid)) ||
        ci->dark_matter.dx_max_sort_old > space_maxreldx * ci->dmin)
        error("Interacting unsorted cells.");
    if (!(cj->dark_matter.sorted & (1 << sid)) ||
        cj->dark_matter.dx_max_sort_old > space_maxreldx * cj->dmin)
        error("Interacting unsorted cells.");

#ifdef SWIFT_DEBUG_CHECKS
    /* Pick-out the sorted lists. */
    const struct sort_entry *restrict sort_i = cell_get_dark_matter_sorts(ci, sid);
    const struct sort_entry *restrict sort_j = cell_get_dark_matter_sorts(cj, sid);

    /* Check that the dx_max_sort values in the cell are indeed an upper
     bound on particle movement. */
      for (int pid = 0; pid < ci->dark_matter.count; pid++) {
        const struct dmpart *p = &ci->dark_matter.parts[sort_i[pid].i];
        if (dmpart_is_inhibited(p, e)) continue;

        const float d = p->x[0] * runner_shift[sid][0] +
                        p->x[1] * runner_shift[sid][1] +
                        p->x[2] * runner_shift[sid][2];
        if (fabsf(d - sort_i[pid].d) - ci->dark_matter.dx_max_sort >
                1.0e-4 * max(fabsf(d), ci->dark_matter.dx_max_sort_old) &&
            fabsf(d - sort_i[pid].d) - ci->dark_matter.dx_max_sort >
                ci->width[0] * 1.0e-10)
          error(
              "particle shift diff exceeds dx_max_sort in cell ci. ci->nodeID=%d "
              "cj->nodeID=%d d=%e sort_i[pid].d=%e ci->dark_matter.dx_max_sort=%e "
              "ci->dark_matter.dx_max_sort_old=%e",
              ci->nodeID, cj->nodeID, d, sort_i[pid].d, ci->dark_matter.dx_max_sort,
              ci->dark_matter.dx_max_sort_old);
      }
      for (int pjd = 0; pjd < cj->dark_matter.count; pjd++) {
        const struct dmpart *p = &cj->dark_matter.parts[sort_j[pjd].i];
        if (dmpart_is_inhibited(p, e)) continue;

        const float d = p->x[0] * runner_shift[sid][0] +
                        p->x[1] * runner_shift[sid][1] +
                        p->x[2] * runner_shift[sid][2];
        if (fabsf(d - sort_j[pjd].d) - cj->dark_matter.dx_max_sort >
                1.0e-4 * max(fabsf(d), cj->dark_matter.dx_max_sort_old) &&
            fabsf(d - sort_j[pjd].d) - cj->dark_matter.dx_max_sort >
                cj->width[0] * 1.0e-10)
          error(
              "particle shift diff exceeds dx_max_sort in cell cj. cj->nodeID=%d "
              "ci->nodeID=%d d=%e sort_j[pjd].d=%e cj->dark_matter.dx_max_sort=%e "
              "cj->dark_matter.dx_max_sort_old=%e",
              cj->nodeID, ci->nodeID, d, sort_j[pjd].d, cj->dark_matter.dx_max_sort,
              cj->dark_matter.dx_max_sort_old);
      }
#endif /* SWIFT_DEBUG_CHECKS */

#ifdef SWIFT_USE_NAIVE_INTERACTIONS
    DOPAIR2_NAIVE(r, ci, cj);
#else
    DOPAIR2(r, ci, cj, sid, shift);
#endif
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
    const struct unit_system *us = e->internal_units;
    const struct sidm_props *sidm_props = e->sidm_properties;
    const int with_cosmology = e->policy & engine_policy_cosmology;
    const double time_base = e->time_base;
    const integertime_t t_current = e->ti_current;

    TIMER_TIC;
    
    struct dmpart *restrict dmparts = c->dark_matter.parts;
    const int count = c->dark_matter.count;
    
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

                /* Compute the pairwise distance. */
                float r2 = 0.0f;
                float dx[3];
                for (int k = 0; k < 3; k++) {
                    dx[k] = pj->x[k] - pix[k];
                    r2 += dx[k] * dx[k];
                }

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
                /* Check that particles have been drifted to the current time */
                if (pi->ti_drift != e->ti_current)
                  error("Particle pi not drifted to current time");
                if (pj->ti_drift != e->ti_current)
                  error("Particle pj not drifted to current time");
#endif

                /* Hit or miss? */
                if (r2 < hig2 || r2 < hj * hj * dm_kernel_gamma2) {
                    
                    runner_iact_nonsym_dark_matter_sidm(r2, dx, hj, hi, pj, pi, a, H, time_base,
                                                        t_current, cosmo, with_cosmology, sidm_props, us);
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

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
                /* Check that particles have been drifted to the current time */
                if (pi->ti_drift != e->ti_current)
                  error("Particle pi not drifted to current time");
                if (pj->ti_drift != e->ti_current)
                  error("Particle pj not drifted to current time");
#endif

                /* Hit or miss? */
                if (r2 < hig2 || r2 < hj * hj * dm_kernel_gamma2) {
                    
                    /* Does pj need to be updated too? */
                    if (dmpart_is_active(pj, e)) {
                        
                        runner_iact_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                                     t_current, cosmo, with_cosmology, sidm_props, us);
                    } else {
                        
                        runner_iact_nonsym_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                                            t_current, cosmo, with_cosmology, sidm_props, us);

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
    const struct unit_system *us = e->internal_units;
    const struct sidm_props *sidm_props = e->sidm_properties;
    const int with_cosmology = e->policy & engine_policy_cosmology;
    const double time_base = e->time_base;
    const integertime_t t_current = e->ti_current;

    TIMER_TIC;

    /* Anything to do here? */
    if (!cell_is_active_dark_matter(c, e)) return;

    /* Cosmological terms */
    const float a = cosmo->a;
    const float H = cosmo->H;
    
    const int count = c->dark_matter.count;
    struct dmpart *restrict dmparts = c->dark_matter.parts;

    /* Loop over the dmparts in ci. */
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

        /* Loop over the dmparts in cj. */
        for (int pjd = pid + 1; pjd < count; pjd++) {
            
            /* Get a pointer to the jth particle. */
            struct dmpart *restrict pj = &dmparts[pjd];
            
            const float hj = pj->h;
            const float hjg2 = hj * hj * dm_kernel_gamma2;
            const int pj_active = dmpart_is_active(pj, e);

            /* Compute the pairwise distance. */
            const float pjx[3] = {(float)(pj->x[0] - c->loc[0]),
                                  (float)(pj->x[1] - c->loc[1]),
                                  (float)(pj->x[2] - c->loc[2])};
            float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
            const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

            const int doi = pi_active && ((r2 < hig2) || (r2 < hjg2));
            const int doj = pj_active && ((r2 < hig2) || (r2 < hjg2));

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
            /* Check that particles have been drifted to the current time */
            if (pi->ti_drift != e->ti_current)
                error("Particle pi not drifted to current time");
            if (pj->ti_drift != e->ti_current)
                error("Particle pj not drifted to current time");
#endif

            /* Hit or miss? */
            if (doi && doj) {
                
                runner_iact_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                             t_current, cosmo, with_cosmology, sidm_props, us);

            } else if (doi) {

                runner_iact_nonsym_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, time_base,
                                                    t_current, cosmo, with_cosmology, sidm_props, us);

            } else if (doj) {
                
                dx[0] = -dx[0];
                dx[1] = -dx[1];
                dx[2] = -dx[2];

                runner_iact_nonsym_dark_matter_sidm(r2, dx, hj, hi, pj, pi, a, H, time_base,
                                                    t_current, cosmo, with_cosmology, sidm_props, us);

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
void DOSELF2_BRANCH(struct runner *r, struct cell *c) {

    const struct engine *restrict e = r->e;

    /* Anything to do here? */
    if (c->dark_matter.count == 0) return;

    /* Anything to do here? */
    if (!cell_is_active_dark_matter(c, e)) return;

    /* Did we mess up the recursion? */
    if (c->dark_matter.h_max_old * dm_kernel_gamma > c->dmin)
        error("Cell smaller than smoothing length");

    /* Check that cells are drifted. */
    if (!cell_are_dmpart_drifted(c, e)) error("Interacting undrifted cell.");

#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
    DOSELF2_NAIVE(r, c);
#else
    DOSELF2(r, c);
#endif
}


/**
 * @brief Compute grouped sub-cell interactions for pairs (symmetric case)
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param cj The second #cell.
 * @param gettimer Do we have a timer ?
 *
 * @todo Hard-code the sid on the recursive calls to avoid the
 * redundant computations to find the sid on-the-fly.
 */
void DOSUB_PAIR2(struct runner *r, struct cell *ci, struct cell *cj,
                 int gettimer) {

    const struct engine *e = r->e;
    struct space *s = e->s;

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
                DOSUB_PAIR2(r, ci->progeny[pid], cj->progeny[pjd], 0);
        }
    }

        /* Otherwise, compute the pair directly. */
    else if (cell_is_active_dark_matter(ci, e) || cell_is_active_dark_matter(cj, e)) {

        /* Make sure both cells are drifted to the current timestep. */
        if (!cell_are_dmpart_drifted(ci, e) || !cell_are_dmpart_drifted(cj, e))
            error("Interacting undrifted cells.");

        /* Do any of the cells need to be sorted first? */
        if (!(ci->dark_matter.sorted & (1 << sid)) ||
            ci->dark_matter.dx_max_sort_old > ci->dmin * space_maxreldx)
            error(
                "Interacting unsorted cell. ci->dark_matter.dx_max_sort_old=%e ci->dmin=%e "
                "ci->sorted=%d sid=%d",
                ci->dark_matter.dx_max_sort_old, ci->dmin, ci->dark_matter.sorted, sid);
        if (!(cj->dark_matter.sorted & (1 << sid)) ||
            cj->dark_matter.dx_max_sort_old > cj->dmin * space_maxreldx)
            error(
                "Interacting unsorted cell. cj->dark_matter.dx_max_sort_old=%e cj->dmin=%e "
                "cj->sorted=%d sid=%d",
                cj->dark_matter.dx_max_sort_old, cj->dmin, cj->dark_matter.sorted, sid);

        /* Compute the interactions. */
        DOPAIR2_BRANCH(r, ci, cj);
    }
    if (gettimer) TIMER_TOC(TIMER_DOSUB_PAIR);
}

/**
 * @brief Compute grouped sub-cell interactions for self tasks (symmetric case)
 *
 * @param r The #runner.
 * @param ci The first #cell.
 * @param gettimer Do we have a timer ?
 */
void DOSUB_SELF2(struct runner *r, struct cell *ci, int gettimer) {

    TIMER_TIC;

    /* Should we even bother? */
    if (ci->dark_matter.count == 0 || !cell_is_active_dark_matter(ci, r->e)) return;

    /* Recurse? */
    if (cell_can_recurse_in_self_dark_matter_task(ci)) {

        /* Loop over all progeny. */
        for (int k = 0; k < 8; k++)
            if (ci->progeny[k] != NULL) {
                DOSUB_SELF2(r, ci->progeny[k], 0);
                for (int j = k + 1; j < 8; j++)
                    if (ci->progeny[j] != NULL)
                        DOSUB_PAIR2(r, ci->progeny[k], ci->progeny[j], 0);
            }

    }
        /* Otherwise, compute self-interaction. */
    else {
        DOSELF2_BRANCH(r, ci);
    }
    if (gettimer) TIMER_TOC(TIMER_DOSUB_SELF);
}



/***********************************************************************/

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
            if (pi->ti_drift != e->ti_current)
                error("Particle pi not drifted to current time");
            if (pj->ti_drift != e->ti_current)
                error("Particle pj not drifted to current time");
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

   /* Get the sorting index. */
   int sid = 0;
   for (int k = 0; k < 3; k++)
      sid = 3 * sid + ((cj->loc[k] - ci->loc[k] + shift[k] < 0) ? 0
                      : (cj->loc[k] - ci->loc[k] + shift[k] > 0) ? 2
                                                                 : 1);

   /* Switch the cells around? */
   /*const int flipped = runner_flip[sid];*/
   sid = sortlistID[sid];

    DOPAIR_SUBSET_NAIVE(r, ci, dmparts_i, ind, count, cj, shift);
}


void DOSUB_SUBSET(struct runner *r, struct cell *ci, struct dmpart *dmparts,
                  int *ind, int count, struct cell *cj, int gettimer) {

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
            DOSUB_SUBSET(r, sub, dmparts, ind, count, NULL, 0);
            for (int j = 0; j < 8; j++)
                if (ci->progeny[j] != sub && ci->progeny[j] != NULL)
                    DOSUB_SUBSET(r, sub, dmparts, ind, count, ci->progeny[j], 0);
        }

            /* Otherwise, compute self-interaction. */
        else
            DOSELF_SUBSET_BRANCH(r, ci, dmparts, ind, count);
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
                    DOSUB_SUBSET(r, ci->progeny[pid], dmparts, ind, count, cj->progeny[pjd], 0);
                if (ci->progeny[pid] != NULL && cj->progeny[pjd] == sub)
                    DOSUB_SUBSET(r, cj->progeny[pjd], dmparts, ind, count, ci->progeny[pid], 0);
            }
        }

            /* Otherwise, compute the pair directly. */
        else if (cell_is_active_dark_matter(ci, e) || cell_is_active_dark_matter(cj, e)) {

            /* Do any of the cells need to be drifted first? */
            if (!cell_are_dmpart_drifted(cj, e)) error("Cell should be drifted!");

            DOPAIR_SUBSET_BRANCH(r, ci, dmparts, ind, cj, count);
        }

    } /* otherwise, pair interaction. */

    if (gettimer) TIMER_TOC(timer_dosub_subset);
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

#ifdef SWIFT_DEBUG_CHECKS
        if (!dmpart_is_active(pi, e)) error("Inactive dmparticle in subset function!");
#endif

        /* Loop over the parts in cj. */
        for (int pjd = 0; pjd < count_i; pjd++) {

            /* Get a pointer to the jth particle. */
            struct dmpart *restrict pj = &dmparts_j[pjd];

            /* Skip oneself */
            if (pi == pj) continue;

            /* Skip inhibited particles. */
            if (dmpart_is_inhibited(pj, e)) continue;

            const float hj = pj->h;

            /* Compute the pairwise distance. */
            const float pjx[3] = {(float)(pj->x[0] - ci->loc[0]),
                                  (float)(pj->x[1] - ci->loc[1]),
                                  (float)(pj->x[2] - ci->loc[2])};
            float dx[3] = {pix[0] - pjx[0], pix[1] - pjx[1], pix[2] - pjx[2]};
            const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

#if defined(SWIFT_DEBUG_CHECKS) && defined(DO_DRIFT_DEBUG_CHECKS)
            /* Check that particles have been drifted to the current time */
            if (pi->ti_drift != e->ti_current)
                error("Particle pi not drifted to current time");
            if (pj->ti_drift != e->ti_current)
                error("Particle pj not drifted to current time");
#endif
            /* Hit or miss? */
            if (r2 < hig2) {

                runner_iact_nonsym_dark_matter_density(r2, dx, hi, hj, pi, pj, a, H);
            }
        } /* loop over the parts in cj. */
    }     /* loop over the parts in ci. */

    TIMER_TOC(timer_doself_subset);
}

/**
 * @brief Determine which version of DOSELF_SUBSET needs to be called depending
 * on the optimisation level.

 * @param r The #runner.
 * @param ci The first #cell.
 * @param parts The #part to interact.
 * @param ind The list of indices of particles in @c ci to interact with.
 * @param count The number of particles in @c ind.
 */
void DOSELF_SUBSET_BRANCH(struct runner *r, struct cell *restrict ci,
                          struct dmpart *restrict dmparts, int *restrict ind,
                          int count) {

    DOSELF_SUBSET(r, ci, dmparts, ind, count);
}



