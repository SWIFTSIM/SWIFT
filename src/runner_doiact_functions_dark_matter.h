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
        
#ifdef WITH_MPI
    struct space *s = e->s;
    struct dmpart *dmparts_foreign = s->dmparts_foreign;
    const size_t nr_dmparts_foreign = s->nr_dmparts_foreign;
#endif

    
    TIMER_TIC;
    
    /* Anything to do here? */
    if (cj->dark_matter.count == 0 || ci->dark_matter.count == 0) return;
    if (!cell_is_active_dark_matter(ci, e)) return;
    
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
        const float hi = pi->sidm_data.h_sidm;
        const float hig2 = hi * hi;
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
            
            const int pj_active = dmpart_is_active(pj, e);
            const float hj = pj->sidm_data.h_sidm;
            const float hjg2 = hj * hj;
            
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
                
                runner_iact_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, dti, dtj, ti_begin, sidm_props, us);
                
            } else if (doi) {
                
                runner_iact_nonsym_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, dti, dtj, ti_begin, sidm_props, us);
                
            } else if (doj) {
                
                dx[0] = -dx[0];
                dx[1] = -dx[1];
                dx[2] = -dx[2];
                
                runner_iact_nonsym_dark_matter_sidm(r2, dx, hj, hi, pj, pi, a, H, dtj, dti, ti_begin, sidm_props, us);
                
            }
            
#ifdef WITH_MPI

            /* If hit then calculate collisions/change of direction and update velocities */
            if (pi->sidm_data.sidm_flag == 1){

                /* First check whether one of the local particles has been hit by foreign part */
                for (size_t pidf = 0; pidf < nr_dmparts_foreign; ++pidf) {
                    
                    struct dmpart *dmpf = &dmparts_foreign[pidf];
                    
                    /* Get a handle on the dm foreign part. */
                    if ((dmpf->id_or_neg_offset == pi->id_or_neg_offset) || (dmpf->id_or_neg_offset == pj->id_or_neg_offset)) {
                        
                        if (dmpf->sidm_data.sidm_flag == 1) {
                        
                        message("Local dmpart has been hit by a foreign one)");
                        
                        }
                    }
                }
            }
#endif

            
        } /* loop over the parts in cj. */
    }   /* loop over the parts in ci. */
        
    TIMER_TOC(TIMER_DOPAIR);
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
        const float hi = pi->sidm_data.h_sidm;
        const float hig2 = hi * hi;
        const float pix[3] = {(float)(pi->x[0] - c->loc[0]),
            (float)(pi->x[1] - c->loc[1]),
            (float)(pi->x[2] - c->loc[2])};
        
        /* Loop over the parts in cj. */
        for (int pjd = pid + 1; pjd < count; pjd++) {
            
            /* Get a pointer to the jth particle. */
            struct dmpart *restrict pj = &dmparts[pjd];
            
            /* Skip inhibited particles. */
            if (dmpart_is_inhibited(pj, e)) continue;
            
            const float hj = pj->sidm_data.h_sidm;
            const float hjg2 = hj * hj;
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
            
            const int doi = pi_active && ((r2 < hig2) || (r2 < hjg2));
            const int doj = pj_active && ((r2 < hig2) || (r2 < hjg2));
            
#ifdef SWIFT_DEBUG_CHECKS
            /* Check that particles have been drifted to the current time */
            if (pi->ti_drift != e->ti_current)
            error("Particle pi not drifted to current time");
            if (pj->ti_drift != e->ti_current)
            error("Particle pj not drifted to current time");
#endif
            
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
        const float hi = pi->sidm_data.h_sidm;
        const float hig2 = hi * hi;
        
        /* Is the ith particle not active? */
        if (!dmpart_is_active(pi, e)) {
            
            /* Loop over the other particles */
            for (int pjd = firstdt; pjd < countdt; pjd++) {
                
                /* Get a pointer to the jth particle. */
                struct dmpart *restrict pj = &dmparts[indt[pjd]];
                const float hj = pj->sidm_data.h_sidm;
                
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
                
#ifdef SWIFT_DEBUG_CHECKS
                /* Check that particles have been drifted to the current time */
                if (pi->ti_drift != e->ti_current)
                error("Particle pi not drifted to current time");
                if (pj->ti_drift != e->ti_current)
                error("Particle pj not drifted to current time");
#endif
                
                /* Hit or miss? */
                if (r2 < hig2 || r2 < hj * hj) {
                    
                    runner_iact_nonsym_dark_matter_sidm(r2, dx, hj, hi, pj, pi, a, H, dtj, dti, ti_begin, sidm_props, us);
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
                
                const float hj = pj->sidm_data.h_sidm;
                
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
                
#ifdef SWIFT_DEBUG_CHECKS
                /* Check that particles have been drifted to the current time */
                if (pi->ti_drift != e->ti_current)
                error("Particle pi not drifted to current time");
                if (pj->ti_drift != e->ti_current)
                error("Particle pj not drifted to current time");
#endif
                
                /* Hit or miss? */
                if (r2 < hig2 || r2 < hj * hj) {
                    
                    /* Does pj need to be updated too? */
                    if (dmpart_is_active(pj, e)) {
                        
                        runner_iact_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, dti, dtj, ti_begin, sidm_props, us);

                    } else {
                        
                        runner_iact_nonsym_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, dti, dtj, ti_begin, sidm_props, us);
                    }
                }
            } /* loop over all other particles. */
        }
    } /* loop over all particles. */
    
    free(indt);
    
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
    
    /* Check that cells are drifted. */
    if (!cell_are_dmpart_drifted(c, e)) error("Interacting undrifted cell.");
    
#if defined(SWIFT_USE_NAIVE_INTERACTIONS)
    DOSELF2_NAIVE(r, c);
#else
    DOSELF2(r, c);
#endif
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
        const float hi = pi->sidm_data.h_sidm;
        const float hig2 = hi * hi;
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
            const float hj = pj->sidm_data.h_sidm;
            const float hjg2 = hj * hj;
            
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
            
#ifdef SWIFT_DEBUG_CHECKS
            /* Check that particles have been drifted to the current time */
            if (pi->ti_drift != e->ti_current)
            error("Particle pi not drifted to current time");
            if (pj->ti_drift != e->ti_current)
            error("Particle pj not drifted to current time");
#endif
            
            /* Hit or miss? */
            if (r2 < hig2 || r2 < hjg2) {
                
                if (pi_active && pj_active) {
                    
                    runner_iact_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, dti, dtj, ti_begin, sidm_props, us);

                } else if (pi_active) {
                    
                    runner_iact_nonsym_dark_matter_sidm(r2, dx, hi, hj, pi, pj, a, H, dti, dtj, ti_begin, sidm_props, us);

                } else if (pj_active) {
                    
                    dx[0] = -dx[0];
                    dx[1] = -dx[1];
                    dx[2] = -dx[2];
                    
                    runner_iact_nonsym_dark_matter_sidm(r2, dx, hj, hi, pj, pi, a, H, dtj, dti, ti_begin, sidm_props, us);
                }
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
void DOPAIR2_BRANCH(struct runner *r, struct cell *ci, struct cell *cj) {
    
    const struct engine *restrict e = r->e;
    
    /* Anything to do here? */
    if (ci->dark_matter.count == 0 || cj->dark_matter.count == 0) return;
    
    /* Anything to do here? */
    if (!cell_is_active_dark_matter(ci, e) && !cell_is_active_dark_matter(cj, e)) return;
    
    /* Check that cells are drifted. */
    if (!cell_are_dmpart_drifted(ci, e) || !cell_are_dmpart_drifted(cj, e))
    error("Interacting undrifted cells.");
    
    
#ifdef SWIFT_USE_NAIVE_INTERACTIONS
    DOPAIR2_NAIVE(r, ci, cj);
#else
    DOPAIR2(r, ci, cj);
#endif
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
void DOSUB_PAIR2(struct runner *r, struct cell *ci, struct cell *cj) {
    
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
    if (ci->split && cj->split) {
        struct cell_split_pair *csp = &cell_split_pairs[sid];
        for (int k = 0; k < csp->count; k++) {
            const int pid = csp->pairs[k].pid;
            const int pjd = csp->pairs[k].pjd;
            if (ci->progeny[pid] != NULL && cj->progeny[pjd] != NULL)
            DOSUB_PAIR2(r, ci->progeny[pid], cj->progeny[pjd]);
        }
    }
    
    /* Otherwise, compute the pair directly. */
    else if (cell_is_active_dark_matter(ci, e) || cell_is_active_dark_matter(cj, e)) {
        
        /* Make sure both cells are drifted to the current timestep. */
        if (!cell_are_dmpart_drifted(ci, e) || !cell_are_dmpart_drifted(cj, e))
        error("Interacting undrifted cells.");
        
        /* Compute the interactions. */
        DOPAIR2_BRANCH(r, ci, cj);
    }
    
    TIMER_TOC(TIMER_DOSUB_PAIR);
}
    
/**
 * @brief Compute grouped sub-cell interactions for self tasks (symmetric case)
 *
 * @param r The #runner.
 * @param ci The first #cell.
 */
void DOSUB_SELF2(struct runner *r, struct cell *ci) {
    
    TIMER_TIC;
    
    /* Should we even bother? */
    if (ci->dark_matter.count == 0 || !cell_is_active_dark_matter(ci, r->e)) return;
    
    /* Recurse? */
    if (ci->split) {
        
        /* Loop over all progeny. */
        for (int k = 0; k < 8; k++)
        if (ci->progeny[k] != NULL) {
            DOSUB_SELF2(r, ci->progeny[k]);
            for (int j = k + 1; j < 8; j++)
            if (ci->progeny[j] != NULL)
            DOSUB_PAIR2(r, ci->progeny[k], ci->progeny[j]);
        }
        
    }
    
    /* Otherwise, compute self-interaction. */
    else {
        DOSELF2_BRANCH(r, ci);
    }
    
    TIMER_TOC(TIMER_DOSUB_SELF);
}

