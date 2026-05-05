/*******************************************************************************
 * This file is part of SWIFT.
 * Copyright (c) 2026 Darwin Roduit (darwin.roduit@epfl.ch)
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
#include <config.h>

/* Some standard headers. */
#include <float.h>
#include <limits.h>
#include <stdlib.h>

/* MPI headers. */
#ifdef WITH_MPI
#include <mpi.h>
#endif

/* This object's header. */
#include "runner.h"

/* Local headers. */
#include "active.h"
#include "cell.h"
#include "engine.h"
#include "error.h"
#include "feedback.h"
#include "fof.h"
#include "gravity.h"
#include "hydro.h"
#include "runner_doiact_sinks.h"
#include "space.h"
#include "star_formation.h"
#include "stars.h"
#include "timers.h"
#include "timestep_limiter.h"
#include "tracers.h"

/**
 * @brief Top-level function for HII ionization feedback.
 *
 * This function recurses down the cell hierarchy until it reaches a level
 * where the cell size is comparable to the maximum HII radius of the stars
 * within it. At that point, it calls the 'branch' function to handle the
 * actual neighbor gathering.
 *
 * @param r The #runner thread.
 * @param c The #cell to process.
 * @param timer If true, records the timing of this operation.
 */
void runner_do_stars_hii_ionization_feedback(struct runner *r, struct cell *c,
                                             int timer) {
#ifdef IONIZATION_FEEDBACK_LOOP
  struct engine *e = r->e;

  /* TODO: Or h_max_old? */
  /* TODO: Which cell h_hii_max shall I look at? The current code changes the
     value as we do deeper in the cell hierarchy. Ideally, we want to use the
     smallest hmax for a cell, process the stars and parts at this level for
     this cell, and go higher for cells/stars with higher h_hii */

  float r_hii_max = c->hydro.super->stars.h_hii_max_old * kernel_gamma;
  if (c->stars.h_hii_max_old <= 0.0) {
    r_hii_max = c->hydro.super->stars.h_max_old * kernel_gamma;
  }

  const float interaction_limit = 1.2f * r_hii_max;
  const int can_recurse = c->split && (interaction_limit < 0.5f * c->dmin);

  if (c->stars.count == 0 || c->hydro.count == 0 || !cell_is_active_stars(c, e))
    return;

#ifdef SWIFT_DEBUG_CHECKS
	if (c->cellID != c->hydro.super->cellID && timer == 1) {
	  warning("Not running the HII reionization task on the hydro super (c = %lld, c->hydro.super = %lld)!", c->cellID, c->hydro.super->cellID);
	}
#endif

  TIMER_TIC;

  /*
   * Decision: Should we process this cell or go deeper?
   * We stop recursing if:
   * 1. The cell is a leaf (!c->split).
   * 2. The cell is small enough that we can't efficiently split the
   *    interaction radius further (similar to Mosaics logic).
   * 3. The R_HII_max covers a significant fraction of the cell.
   */
  /* Is the cell split and not smaller than the smoothing length? */
  if (can_recurse) {
    /* Keep recursing deeper into the super-cell hierarchy. */
    for (int k = 0; k < 8; k++) {
      if (c->progeny[k] != NULL) {
	struct cell *restrict cp = c->progeny[k];
        runner_do_stars_hii_ionization_feedback(r, cp, 0);
	c->stars.h_hii_max = max(c->stars.h_hii_max, cp->stars.h_hii_max);
      }
    }
  } else {
    /* We have reached the 'Working Level' */
    runner_do_stars_hii_ionization_feedback_branch(r, c);
  }

  if (timer) TIMER_TOC(timer_stars_hii_ionization_feedback);
#endif
}

/**
 * @brief Branching point to coordinate HII ionization for stars in a cell.
 *
 * This function manages the neighbor gathering for all active stars in the
 * cell. It allocates a temporary buffer, gathers neighbors from the cell itself
 * (self) and its neighbors (pair) by climbing the cell hierarchy, and
 * prepares the gathered particles for sorting and ionization.
 *
 * @param r The #runner thread.
 * @param c The #cell containing the stars to process.
 */
void runner_do_stars_hii_ionization_feedback_branch(struct runner *r,
                                                    struct cell *c) {

  struct engine *e = r->e;
  const struct cosmology *cosmo = e->cosmology;
  /* const struct feedback_props *fb_props = e->feedback_props; */
  const struct hydro_props *hydro_props = e->hydro_properties;
  const struct phys_const *phys_const = e->physical_constants;
  const struct unit_system *us = e->internal_units;
  const struct cooling_function_data *cooling = e->cooling_func;

  struct spart *restrict sparts = c->stars.parts;
  const int scount = c->stars.count;

  float r_hii_max = c->hydro.super->stars.h_hii_max_old * kernel_gamma;
  if (c->stars.h_hii_max_old <= 0.0) {
    r_hii_max = c->hydro.super->stars.h_max_old * kernel_gamma;
  }

  const float interaction_limit = 1.2f * r_hii_max;

  if (c->stars.count == 0 || c->hydro.count == 0 || !cell_is_active_stars(c, e))
    return;

#ifdef SWIFT_DEBUG_CHECKS
  /* Did we mess up the recursion? */
  if (interaction_limit > c->dmin)
    error("Cell smaller than HII interaction length");
#endif

  /* message("[c = %lld, super = %lld] interaction_limit = %e, cell_dmin = %e",
   */
  /*         c->cellID, c->hydro.super->cellID, interaction_limit, c->dmin); */

  /* TODO: Handle the case where this value is too small. */
  /* TODO: Add multiple tries if the star has not exhausted its photons */
  const int max_ngbs = 128;
  struct hii_neighbor *ngb_buffer =
      malloc(max_ngbs * sizeof(struct hii_neighbor));

  for (int sid = 0; sid < scount; sid++) {

    /* Get a hold of the ith spart in ci. */
    struct spart *si = &sparts[sid];

    /* Is this part within the timestep? */
    if (!spart_is_active(si, e) && !feedback_is_active(si, e) && !feedback_is_HII_ionization_active(si, e))
      continue;

#ifdef SWIFT_DEBUG_CHECKS
    /* Check that particles have been drifted to the current time */
    if (si->ti_drift != e->ti_current)
      error("Particle si (%lld) not drifted to current time c = %lld, "
              "c->super = %lld", si->id, c->cellID, c->super->cellID);

#endif

    /* Logic: If the buffer was full, there might be more neighbors just 
     * outside the current R_max of the buffer. We retry until the buffer 
     * is no longer maxed out OR the star runs out of photons. */
    int buffer_was_full = 0;
    do {
      int count_found = 0;
      buffer_was_full = 0;

      /***************************************************/
      /* First loop over particles in the current cell */
      runner_do_stars_hii_ionization_feedback_self(r, c, si, ngb_buffer, max_ngbs,
						   &count_found);

      /***************************************************/
      /* Now loop over particles in the neighboring cells */

      /* Climb up the cell hierarchy. */
      for (struct cell *finger = c; finger != NULL; finger = finger->parent) {
	/* These are defined at the super level... When we reach the progeny,
	   these will be NULL... So we need to grab the finger first */
	for (struct link *l = finger->stars.density; l != NULL; l = l->next) {
	  /* We have already handled the self case */
	  if (l->t->type == task_type_self) continue;

	  struct cell *c_in = (l->t->cj->hydro.super == c->hydro.super) ? l->t->ci->hydro.super : l->t->cj->hydro.super;

#ifdef SWIFT_DEBUG_CHECKS
	  if (c_in->hydro.super->cellID == c->hydro.super->cellID) {
	    warning("cj (%lld) has the same hydro super (%lld) than me (%lld)!", c_in->cellID, c->hydro.super->cellID, c->cellID);
	  }
#endif

	  struct cell *cj = l->t->cj;
	  struct cell *ci = l->t->ci;
	  message("[%lld, %lld], hydro super = %lld , %lld", ci->cellID, cj->cellID, ci->hydro.super->cellID, cj->hydro.super->cellID);

	  /* Now, find the correct level... */
	  const int can_recurse_j =
            c_in->split && (interaction_limit < 0.5f * c_in->dmin);
	  if (can_recurse_j) {
	    /* Keep recursing deeper into the super-cell hierarchy. */
	    for (int k = 0; k < 8; k++)
	      if (c_in->progeny[k] != NULL)
		runner_do_stars_hii_ionization_feedback_pair(r, c, c_in->progeny[k], si, ngb_buffer, max_ngbs,
							     &count_found);
	  } else {
	    /* We have reached the 'Working Level' */
	    runner_do_stars_hii_ionization_feedback_pair(r, c, c_in, si, ngb_buffer, max_ngbs, &count_found);
	  } /* Recurse */
	} /* Neighbour search */
      } /* Climb up in the cell hierarchy */

      /* Flag if we hit the limit */
      if (count_found == max_ngbs) buffer_was_full = 1;    

      /***************************************************/
      /* It's time to sort the gas particles */
      if (count_found > 0) {
#ifdef SWIFT_DEBUG_CHECKS
	/* Verify that the neighbors are properly sorted by distance */
	for (int k = 0; k < count_found - 1; k++) {
	  if (ngb_buffer[k].r2 > ngb_buffer[k + 1].r2) {
	    error(
		  "HII neighbor buffer not properly sorted! "
		  "Index %d (r2=%e) is larger than index %d (r2=%e).",
		  k, ngb_buffer[k].r2, k + 1, ngb_buffer[k + 1].r2);
	  }
	}
#endif
	const integertime_t ti_begin = get_integer_time_begin(e->ti_current - 1, si->time_bin);

	/* Now let's ionize the gas particles! */
	for (int k = 0; k < count_found; k++) {

	  /* No more photons to consume */
	  if (si->feedback_data.radiation.dot_N_ion <= 0) {
	    message("Star has exhausted all its ionizing photons!");
            break;
	  }
        
	  struct part *pj = ngb_buffer[k].p;
	  struct xpart *xpj = ngb_buffer[k].xp;

	  /* Do the ionization */
	  /* 1. Tag gas as ionized (atomics)
	     2. Flag to be synchronized (atomics)
	     3. Consume photons from the star
	     4. Compute M_HII and r_HII for the star
	     5. Update the stars r_HII.
	     6. Update the cell max r_HII after finising processing all star
	     particles
	  */
	  /* feedback_do_HII_ionization(p, xp); */

	  /* Fast non-atomic check: if already ionized, just move on. */
          if (xpj->tracers_data.HII_region.is_ionized) continue;

	  message("Ionize %lld (budget = %e)!", pj->id, radiation_get_star_ionization_rate(si));

	  const double Delta_dot_N_ion = radiation_get_part_rate_to_fully_ionize(phys_const, hydro_props, us, cosmo, cooling, pj, xpj);

	  /* Case 1: Ionization is guaranteed */
	  if (Delta_dot_N_ion <= radiation_get_star_ionization_rate(si)) {
	    if (atomic_cas(&xpj->tracers_data.HII_region.is_ionized, 0, 1) == 0) {

	      /* Flag the particle to be synchronized on the timeline.
		 We still use atomic_or here because other stars might be
		 tripping the limiter in feedback loop. */
	      atomic_or(&pj->limiter_data.to_be_synchronized, 1);

	      /* Add the star ID */
	      xpj->tracers_data.HII_region.star_id = si->id;

	      /* Consume photons from the star */
	      radiation_consume_ionizing_photons(si, Delta_dot_N_ion);

	      /* Update HII region properties */
	      si->feedback_data.radiation.mass_HII_region += hydro_get_mass(pj);
	      si->h_hii = max(si->h_hii, sqrtf(ngb_buffer[k].r2) * kernel_gamma_inv);
	    }
	    /* If CAS failed, someone else grabbed it; just continue. */
	  } else {

	    /* If we cannot fully ionize, compute a probability to determine if
	       we fully ionize pj or not and draw the random number.  */
	    const double dot_N_ion = radiation_get_star_ionization_rate(si);
	    const float proba = dot_N_ion / Delta_dot_N_ion;
	    const float random_number = random_unit_interval(si->id, ti_begin, random_number_HII_regions);

	    /* If we are lucky, do the ionization */
	    if (random_number <= proba) {
	      /* We won the roll! Now try to claim the particle. */
	      if (atomic_cas(&xpj->tracers_data.HII_region.is_ionized, 0, 1) == 0) {
		atomic_or(&pj->limiter_data.to_be_synchronized, 1);

		xpj->tracers_data.HII_region.star_id = si->id;

		/* The star is now empty of photons */
		radiation_consume_ionizing_photons(si, Delta_dot_N_ion);
		si->feedback_data.radiation.mass_HII_region += hydro_get_mass(pj);
		si->h_hii = max(si->h_hii, sqrtf(ngb_buffer[k].r2) * kernel_gamma_inv);

		/* We are out of photons, stop looking for neighbors for this star */
		break;
	      }
	    } else {
	      /* We lost the roll. We still consume the remaining photons. */
	      radiation_consume_ionizing_photons(si, Delta_dot_N_ion);
	      /* Star is empty, stop. */
	      break;
	    }
	  } /* End of probability handling */
	} /* Loop over the sorted particles */
      }
      /* If the star still has photons and we previously filled the buffer, 
       * we must go again to find the neighbors that were 'bumped out'. */
    } while (buffer_was_full && radiation_get_star_ionization_rate(si) > 0);

    c->stars.h_hii_max = max(c->stars.h_hii_max, si->h_hii);
  } /* Loop over sparts */
  free(ngb_buffer);
}

/**
 * @brief Gather gas particles within the HII radius from the star's own cell.
 *
 * Performs a brute-force search over all hydro particles in cell @p c that
 * fall within the HII smoothing length of the star @p si. Candidates are
 * added to the @p buffer if they are not already inhibited or ionized.
 *
 * @param r The #runner thread.
 * @param c The #cell containing the gas particles.
 * @param si The #spart (star) performing the feedback.
 * @param buffer The #hii_neighbor array to store found candidates.
 * @param max_size The maximum capacity of the neighbor buffer.
 * @param count_found (return) The number of neighbors successfully gathered.
 */
void runner_do_stars_hii_ionization_feedback_self(
    struct runner *r, struct cell *c, struct spart *si,
    struct hii_neighbor *buffer, int max_size, int *count_found) {

  struct engine *e = r->e;
  struct part *restrict parts = c->hydro.parts;
  struct xpart *restrict xparts = c->hydro.xparts;
  const int count = c->hydro.count;
  const float six[3] = {(float)(si->x[0] - c->loc[0]),
                        (float)(si->x[1] - c->loc[1]),
                        (float)(si->x[2] - c->loc[2])};

  /* Loop over the parts in c. */
  for (int pjd = 0; pjd < count; pjd++) {

    /* Get a pointer to the jth particle. */
    struct part *restrict pj = &parts[pjd];
    struct xpart *restrict xpj = &xparts[pjd];

    /* const float hj = pj->h; */

    /* Early abort? */
    if (part_is_inhibited(pj, e)) continue;
    if (radiation_is_part_tagged_as_ionized(pj, xpj)) continue;

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (pj->ti_drift != e->ti_current)
          error(
              "Particle pj (%lld) not drifted to current time. c = %lld, "
              "c->super = %lld",
              pj->id, c->cellID, c->super->cellID);
#endif

    /* Compute the pairwise distance. */
    const float pjx[3] = {(float)(pj->x[0] - c->loc[0]),
                          (float)(pj->x[1] - c->loc[1]),
                          (float)(pj->x[2] - c->loc[2])};
    const float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
    const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

    runner_hii_buffer_insert(buffer, max_size, count_found, r2, pj, xpj, c);
    /* TODO: If we reach the max_size, we shoudl print a warning. Maybe add a
       flag to increase the stack size? Or to iterate again later by discarding
       these particles, since they will be tagged as ionized */
  } /* Loop in current cell */
}

/**
 * @brief Gather gas particles within the HII radius from a neighboring cell.
 *
 * Similar to the 'self' version, but handles the distance calculations between
 * two different cells (@p ci and @p cj), including periodic boundary
 * conditions/wrapping.
 *
 * @param r The #runner thread.
 * @param ci The #cell containing the star.
 * @param cj The #cell containing the gas particles.
 * @param si The #spart (star) performing the feedback.
 * @param buffer The #hii_neighbor array to store found candidates.
 * @param max_size The maximum capacity of the neighbor buffer.
 * @param count_found (return) The number of neighbors successfully gathered.
 */
void runner_do_stars_hii_ionization_feedback_pair(
    struct runner *r, struct cell *ci, struct cell *cj, struct spart *si,
    struct hii_neighbor *buffer, int max_size, int *count_found) {

  struct engine *e = r->e;
  struct part *restrict parts_j = cj->hydro.parts;
  struct xpart *restrict xparts_j = cj->hydro.xparts;
  const int count_j = cj->hydro.count;

  /* OR h_max_old? */
  const float r_hii_max = ci->stars.h_hii_max_old * kernel_gamma;
  const float interaction_limit = 1.2f * r_hii_max;

#ifdef SWIFT_DEBUG_CHECKS
  /* Did we mess up the recursion? */
  if (interaction_limit > cj->dmin)
    error("Cell smaller than HII interaction length");

  /* Call the function that will do the work */
  if (interaction_limit > cj->dmin) {
    warning(
        "[c = %lld, cj = %lld, star = %lld] cj size is too small. We need to "
        "go up in the cell hierarchy!",
        ci->cellID, cj->cellID, si->id);
  }

  /* message( */
  /*     "[c = %lld, cj = %lld, cj->super = %lld, star = %lld] Neighbour search,
   * " */
  /*     "interaction_limit = %e, cell_dmin = %e", */
  /*     ci->cellID, cj->cellID, cj->hydro.super->cellID, si->id, */
  /*     interaction_limit, cj->dmin); */
#endif

  /* Get the relative distance between the pairs, wrapping. */
  double shift[3] = {0.0, 0.0, 0.0};
  for (int k = 0; k < 3; k++) {
    if (cj->loc[k] - ci->loc[k] < -e->s->dim[k] / 2)
      shift[k] = e->s->dim[k];
    else if (cj->loc[k] - ci->loc[k] > e->s->dim[k] / 2)
      shift[k] = -e->s->dim[k];
  }

  const float six[3] = {(float)(si->x[0] - (cj->loc[0] + shift[0])),
                        (float)(si->x[1] - (cj->loc[1] + shift[1])),
                        (float)(si->x[2] - (cj->loc[2] + shift[2]))};

  for (int pjd = 0; pjd < count_j; pjd++) {

    /* Get a pointer to the jth particle. */
    struct part *restrict pj = &parts_j[pjd];
    struct xpart *restrict xpj = &xparts_j[pjd];

    /* Early abort? */
    if (part_is_inhibited(pj, e)) continue;
    if (radiation_is_part_tagged_as_ionized(pj, xpj)) continue;

#ifdef SWIFT_DEBUG_CHECKS
        /* Check that particles have been drifted to the current time */
        if (pj->ti_drift != e->ti_current)
          error(
              "Particle pj (%lld) not drifted to current time. c = %lld, "
              "c->super = %lld",
              pj->id, cj->cellID, cj->super->cellID);
#endif

    /* Compute the pairwise distance. */
    const float pjx[3] = {(float)(pj->x[0] - cj->loc[0]),
                          (float)(pj->x[1] - cj->loc[1]),
                          (float)(pj->x[2] - cj->loc[2])};
    const float dx[3] = {six[0] - pjx[0], six[1] - pjx[1], six[2] - pjx[2]};
    const float r2 = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];

    runner_hii_buffer_insert(buffer, max_size, count_found, r2, pj, xpj, cj);
    /* TODO: If we reach the max_size, we shoudl print a warning. Maybe add a
       flag to increase the stack size? Or to iterate again later by discarding
       these particles, since they will be tagged as ionized */
  }
}


/**
 * @brief Maintain a sorted buffer by inserting a new neighbor at the correct position.
 * 
 * If the buffer is full, it replaces the furthest element if the new one is closer.
 */
__attribute__((always_inline)) INLINE void runner_hii_buffer_insert(struct hii_neighbor *buffer, int max_size,
                                     int *count_found, float r2, struct part *p,
                                     struct xpart *xp, struct cell *c) {
  
  /* Case A: Buffer is not yet full */
  if (*count_found < max_size) {
    int i = *count_found - 1;
    /* Shift elements to make room (standard insertion) */
    while (i >= 0 && buffer[i].r2 > r2) {
      buffer[i + 1] = buffer[i];
      i--;
    }
    buffer[i + 1].r2 = r2;
    buffer[i + 1].p = p;
    buffer[i + 1].xp = xp;
#ifdef SWIFT_DEBUG_CHECKS
    buffer[i + 1].c = c;
#endif
    (*count_found)++;
  } 
  /* Case B: Buffer is full, check if new particle is closer than the furthest */
  else if (r2 < buffer[max_size - 1].r2) {
    int i = max_size - 2;
    /* Shift elements to replace the furthest */
    while (i >= 0 && buffer[i].r2 > r2) {
      buffer[i + 1] = buffer[i];
      i--;
    }
    buffer[i + 1].r2 = r2;
    buffer[i + 1].p = p;
    buffer[i + 1].xp = xp;
#ifdef SWIFT_DEBUG_CHECKS
    buffer[i + 1].c = c;
#endif
  }
}
