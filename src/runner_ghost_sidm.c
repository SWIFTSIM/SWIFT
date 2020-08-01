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

/* Config parameters. */
#include "../config.h"

/* This object's header. */
#include "runner.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "space_getsid.h"
#include "timers.h"
#include "timestep_limiter.h"

/**
 * @brief Intermediate task after the density to check that the search radius for
 * self-interacting dark matter is correct.
 *
 * @param r The runner thread.
 * @param c The cell.
 */
void runner_do_dark_matter_density_ghost(struct runner *r, struct cell *c) {
    
    struct gpart *restrict gparts = c->grav.parts;
    const struct engine *e = r->e;
    const struct cosmology *cosmo = e->cosmology;
    const float dark_matter_h_max = e->sidm_properties->h_sidm_max;
    const float dark_matter_h_min = e->sidm_properties->h_sidm_min;
    const float eps = e->sidm_properties->h_tolerance;
    const float dark_matter_eta_dim = pow_dimension(e->sidm_properties->eta_neighbours);
    const int max_smoothing_iter = e->sidm_properties->max_smoothing_iterations;
    int redo = 0, count = 0;
    
    /* Running value of the maximal smoothing length */
    double h_max = c->grav.h_max;
    
    TIMER_TIC;
    
    /* Anything to do here? */
    if (c->black_holes.count == 0) return;
    if (!cell_is_active_black_holes(c, e)) return;
    
    /* Recurse? */
    if (c->split) {
        for (int k = 0; k < 8; k++) {
            if (c->progeny[k] != NULL) {
                runner_do_black_holes_density_ghost(r, c->progeny[k], 0);
                
                /* Update h_max */
                h_max = max(h_max, c->progeny[k]->black_holes.h_max);
            }
        }
    } else {
        
        /* Init the list of active particles that have to be updated. */
        int *sid = NULL;
        float *h_0 = NULL;
        float *left = NULL;
        float *right = NULL;
        if ((sid = (int *)malloc(sizeof(int) * c->black_holes.count)) == NULL)
        error("Can't allocate memory for sid.");
        if ((h_0 = (float *)malloc(sizeof(float) * c->black_holes.count)) == NULL)
        error("Can't allocate memory for h_0.");
        if ((left = (float *)malloc(sizeof(float) * c->black_holes.count)) == NULL)
        error("Can't allocate memory for left.");
        if ((right = (float *)malloc(sizeof(float) * c->black_holes.count)) == NULL)
        error("Can't allocate memory for right.");
        for (int k = 0; k < c->black_holes.count; k++)
        if (bpart_is_active(&bparts[k], e)) {
            sid[bcount] = k;
            h_0[bcount] = bparts[k].h;
            left[bcount] = 0.f;
            right[bcount] = black_holes_h_max;
            ++bcount;
        }
        
        /* While there are particles that need to be updated... */
        for (int num_reruns = 0; bcount > 0 && num_reruns < max_smoothing_iter;
             num_reruns++) {
            
            /* Reset the redo-count. */
            redo = 0;
            
            /* Loop over the remaining active parts in this cell. */
            for (int i = 0; i < bcount; i++) {
                
                /* Get a direct pointer on the part. */
                struct bpart *bp = &bparts[sid[i]];
                
#ifdef SWIFT_DEBUG_CHECKS
                /* Is this part within the timestep? */
                if (!bpart_is_active(bp, e))
                error("Ghost applied to inactive particle");
#endif
                
                /* Get some useful values */
                const float h_init = h_0[i];
                const float h_old = bp->h;
                const float h_old_dim = pow_dimension(h_old);
                const float h_old_dim_minus_one = pow_dimension_minus_one(h_old);
                
                float h_new;
                int has_no_neighbours = 0;
                
                if (bp->density.wcount == 0.f) { /* No neighbours case */
                    
                    /* Flag that there were no neighbours */
                    has_no_neighbours = 1;
                    
                    /* Double h and try again */
                    h_new = 2.f * h_old;
                    
                } else {
                    
                    /* Finish the density calculation */
                    black_holes_end_density(bp, cosmo);
                    
                    /* Compute one step of the Newton-Raphson scheme */
                    const float n_sum = bp->density.wcount * h_old_dim;
                    const float n_target = black_holes_eta_dim;
                    const float f = n_sum - n_target;
                    const float f_prime =
                    bp->density.wcount_dh * h_old_dim +
                    hydro_dimension * bp->density.wcount * h_old_dim_minus_one;
                    
                    /* Improve the bisection bounds */
                    if (n_sum < n_target)
                    left[i] = max(left[i], h_old);
                    else if (n_sum > n_target)
                    right[i] = min(right[i], h_old);
                    
#ifdef SWIFT_DEBUG_CHECKS
                    /* Check the validity of the left and right bounds */
                    if (left[i] > right[i])
                    error("Invalid left (%e) and right (%e)", left[i], right[i]);
#endif
                    
                    /* Skip if h is already h_max and we don't have enough neighbours
                     */
                    /* Same if we are below h_min */
                    if (((bp->h >= black_holes_h_max) && (f < 0.f)) ||
                        ((bp->h <= black_holes_h_min) && (f > 0.f))) {
                        
                        black_holes_reset_feedback(bp);
                        
                        /* Ok, we are done with this particle */
                        continue;
                    }
                    
                    /* Normal case: Use Newton-Raphson to get a better value of h */
                    
                    /* Avoid floating point exception from f_prime = 0 */
                    h_new = h_old - f / (f_prime + FLT_MIN);
                    
                    /* Be verbose about the particles that struggle to converge */
                    if (num_reruns > max_smoothing_iter - 10) {
                        
                        message(
                                "Smoothing length convergence problem: iter=%d p->id=%lld "
                                "h_init=%12.8e h_old=%12.8e h_new=%12.8e f=%f f_prime=%f "
                                "n_sum=%12.8e n_target=%12.8e left=%12.8e right=%12.8e",
                                num_reruns, bp->id, h_init, h_old, h_new, f, f_prime, n_sum,
                                n_target, left[i], right[i]);
                    }
                    
                    /* Safety check: truncate to the range [ h_old/2 , 2h_old ]. */
                    h_new = min(h_new, 2.f * h_old);
                    h_new = max(h_new, 0.5f * h_old);
                    
                    /* Verify that we are actually progrssing towards the answer */
                    h_new = max(h_new, left[i]);
                    h_new = min(h_new, right[i]);
                }
                
                /* Check whether the particle has an inappropriate smoothing length
                 */
                if (fabsf(h_new - h_old) > eps * h_old) {
                    
                    /* Ok, correct then */
                    
                    /* Case where we have been oscillating around the solution */
                    if ((h_new == left[i] && h_old == right[i]) ||
                        (h_old == left[i] && h_new == right[i])) {
                        
                        /* Bissect the remaining interval */
                        bp->h = pow_inv_dimension(
                                                  0.5f * (pow_dimension(left[i]) + pow_dimension(right[i])));
                        
                    } else {
                        
                        /* Normal case */
                        bp->h = h_new;
                    }
                    
                    /* If below the absolute maximum, try again */
                    if (bp->h < black_holes_h_max && bp->h > black_holes_h_min) {
                        
                        /* Flag for another round of fun */
                        sid[redo] = sid[i];
                        h_0[redo] = h_0[i];
                        left[redo] = left[i];
                        right[redo] = right[i];
                        redo += 1;
                        
                        /* Re-initialise everything */
                        black_holes_init_bpart(bp);
                        
                        /* Off we go ! */
                        continue;
                        
                    } else if (bp->h <= black_holes_h_min) {
                        
                        /* Ok, this particle is a lost cause... */
                        bp->h = black_holes_h_min;
                        
                    } else if (bp->h >= black_holes_h_max) {
                        
                        /* Ok, this particle is a lost cause... */
                        bp->h = black_holes_h_max;
                        
                        /* Do some damage control if no neighbours at all were found */
                        if (has_no_neighbours) {
                            black_holes_bpart_has_no_neighbours(bp, cosmo);
                        }
                        
                    } else {
                        error(
                              "Fundamental problem with the smoothing length iteration "
                              "logic.");
                    }
                }
                
                /* We now have a particle whose smoothing length has converged */
                
                black_holes_reset_feedback(bp);
                
                /* Check if h_max has increased */
                h_max = max(h_max, bp->h);
            }
            
            /* We now need to treat the particles whose smoothing length had not
             * converged again */
            
            /* Re-set the counter for the next loop (potentially). */
            bcount = redo;
            if (bcount > 0) {
                
                /* Climb up the cell hierarchy. */
                for (struct cell *finger = c; finger != NULL; finger = finger->parent) {
                    
                    /* Run through this cell's density interactions. */
                    for (struct link *l = finger->black_holes.density; l != NULL;
                         l = l->next) {
                        
#ifdef SWIFT_DEBUG_CHECKS
                        if (l->t->ti_run < r->e->ti_current)
                        error("Density task should have been run.");
#endif
                        
                        /* Self-interaction? */
                        if (l->t->type == task_type_self)
                        runner_doself_subset_branch_bh_density(r, finger, bparts, sid,
                                                               bcount);
                        
                        /* Otherwise, pair interaction? */
                        else if (l->t->type == task_type_pair) {
                            
                            /* Left or right? */
                            if (l->t->ci == finger)
                            runner_dopair_subset_branch_bh_density(r, finger, bparts, sid,
                                                                   bcount, l->t->cj);
                            else
                            runner_dopair_subset_branch_bh_density(r, finger, bparts, sid,
                                                                   bcount, l->t->ci);
                        }
                        
                        /* Otherwise, sub-self interaction? */
                        else if (l->t->type == task_type_sub_self)
                        runner_dosub_subset_bh_density(r, finger, bparts, sid, bcount,
                                                       NULL, 1);
                        
                        /* Otherwise, sub-pair interaction? */
                        else if (l->t->type == task_type_sub_pair) {
                            
                            /* Left or right? */
                            if (l->t->ci == finger)
                            runner_dosub_subset_bh_density(r, finger, bparts, sid, bcount,
                                                           l->t->cj, 1);
                            else
                            runner_dosub_subset_bh_density(r, finger, bparts, sid, bcount,
                                                           l->t->ci, 1);
                        }
                    }
                }
            }
        }
        
        if (bcount) {
            error("Smoothing length failed to converge on %i particles.", bcount);
        }
        
        /* Be clean */
        free(left);
        free(right);
        free(sid);
        free(h_0);
    }
    
    /* Update h_max */
    c->black_holes.h_max = h_max;
    
    /* The ghost may not always be at the top level.
     * Therefore we need to update h_max between the super- and top-levels */
    if (c->black_holes.density_ghost) {
        for (struct cell *tmp = c->parent; tmp != NULL; tmp = tmp->parent) {
            atomic_max_f(&tmp->black_holes.h_max, h_max);
        }
    }
    
    if (timer) TIMER_TOC(timer_do_black_holes_ghost);
}

